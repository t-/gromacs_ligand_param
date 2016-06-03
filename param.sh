#!/bin/bash
# requirements:
# babel installed
# gaussian09 installed and sourced as "g09"

LIGANDLIB=./ligands
FFOUTPATH=./ligandsFFOut
GAUSSIAN_RC_FILE=~/.bashrc_g09
GROMACS_RC_FILE=/usr/local/gromacs/GMXRC51
PATH_TO_ACEPYPE=~/Software/acpype-read-only/acpype.py
BASEDIR=$(pwd)
MINIMALBASISSET_INACCURATE=true

# set GMX ff path
export GMXLIB=/netmount/projects/async/forcefields/
# source gaussian environment... so we can call 'g09'
source $GAUSSIAN_RC_FILE
# source gromacs environment... so we can run editconf
source $GROMACS_RC_FILE

function write_gaussian_input {
  local fname=$1
  # create Gaussian input
  echo '%nprocshared=4' > $fname
  echo '%mem=2000MB' >> $fname
  echo '%chk=mol.chk' >> $fname
  # optionally use very small basis set
  if ! $MINIMALBASISSET_INACCURATE; then
    echo 'using 6-31g(d)'
    echo '#opt hf/6-31g(d) iop(6/33=2,6/41=10,6/42=17) pop=mk scf=tight test' >> $fname
  elif $MINIMALBASISSET_INACCURATE; then
    echo 'using STO-3G (dont do this its terrible but fast)'
    echo '#opt hf/sto-3g iop(6/33=2,6/41=10,6/42=17) pop=mk scf=tight test' >> $fname
  fi
  echo -e '\nTitle Card Required\n' >> $fname
  echo "$MOL_CHARGE 1" >> $fname
  #skip header of .com file (first 6 lines)
  tail -n +6 ${MOLNAME}.com  >> $fname
}

function generate_resgen_input {
  local fname=$1
  touch $fname
  echo "INPUT_FILE mol.ac
CONF_NUM 1
ESP_FILE mol.esp
NET_CHARGE $MOL_CHARGE
PREP_FILE: aa.prep
RESIDUE_FILE_NAME: aa.res
RESIDUE_SYMBOL: AA" > $fname
}

function generate_amber_params {
  antechamber -i mol.out -fi gout -o mol.ac -fo ac -at gaff -nc $MOL_CHARGE -c resp -pf y > log_antechamber1.log 2>&1
  antechamber -i mol.ac -fi ac -o mol.mol2 -fo mol2 -pf y > log_antechamber2.log 2>&1
  espgen -i mol.out -o mol.esp > log_espgen.log 2>&1
  residuegen resgen.DYE > log_residuegen.log 2>&1
  parmchk -i mol.ac -o frcmod -f ac > log_parmchk.log 2>&1
  teLeap -f leap.in > log_teLeap.log 2>&1
}

function modify_amber99sb_ff_for_ligand {
  cp -r $BASEDIR/tools/amber99sb-ildn.ff/* .
  cat bon.template > ffbonded.itp
  cat $BASEDIR/$FFOUTPATH/$MOLNAME/top2itp/merged_bonds.itp | egrep -v '\[|;' >> ffbonded.itp
  cat const.template >> ffbonded.itp
  cat ang.template >> ffbonded.itp
  cat $BASEDIR/$FFOUTPATH/$MOLNAME/top2itp/merged_angles.itp | egrep -v '\[|;' >> ffbonded.itp
  cat impr.template >> ffbonded.itp
  cat $BASEDIR/$FFOUTPATH/$MOLNAME/top2itp/merged_impropers.itp | egrep -v '\[|;'>> ffbonded.itp
  cat dihe.template >> ffbonded.itp
  cat $BASEDIR/$FFOUTPATH/$MOLNAME/top2itp/merged_dihedrals.itp | egrep -v '\[|;' >> ffbonded.itp
  cat def.template >> ffbonded.itp
  cat $BASEDIR/$FFOUTPATH/$MOLNAME/top2itp/merged_nbonds.itp >> ffnonbonded.itp
  cat $BASEDIR/$FFOUTPATH/$MOLNAME/top2itp/res_aminoacids.rtp | sed 's/MOL/AA/g' >> aminoacids.rtp
  cat $BASEDIR/$FFOUTPATH/$MOLNAME/top2itp/merged_atomtypes.atp >> atomtypes.atp
}

for PDBFILENAME in $(ls ligands| grep pdb);
do
  MOLNAME=$(basename -s .pdb $PDBFILENAME)
  MOL_CHARGE=$(echo $MOLNAME|tr '_' '\n'| tail -n 1)
  mkdir -p $BASEDIR/$FFOUTPATH/$MOLNAME
  cd $BASEDIR/$FFOUTPATH/$MOLNAME
  echo -e "\n\n NEW LIGAND \n processing $MOLNAME $(pwd)\n\n"

  #prepare input files and run gaussian optimization
  mkdir -p gau
  cp $BASEDIR/$LIGANDLIB/$MOLNAME.pdb gau
  cd gau
  # convert .pdb to .com file
  babel ${MOLNAME}.pdb ${MOLNAME}.com > log_babel.log 2>&1

  write_gaussian_input mol.com

  echo 'running optimization - this will take some time'
  if [ ! -e mol.out ]; then
    g09 < mol.com > mol.out
  fi

  echo 'finished optimization, checking if job was successful'
  gauoptstate=$(tail -n 1 mol.out | awk '{print $1}')
  if [ ! $gauoptstate == 'Normal' ];then
    echo "gaussian run crashed in $(pwd), continue?"
    select yn in "Yes" "No"; do
        case $yn in
            Yes ) break;;
            No ) exit;;
        esac
    done
  fi

  # parametrize ligand
  mkdir -p $BASEDIR/$FFOUTPATH/$MOLNAME/ligand
  cd $BASEDIR/$FFOUTPATH/$MOLNAME/ligand
  cp $BASEDIR/$FFOUTPATH/$MOLNAME/gau/mol.out .
  cp $BASEDIR/tools/leap.in .

  #generate resgen.DYE file
  generate_resgen_input resgen.DYE

  echo 'trying to auto generate GAFF parameters'
  # trying to do the rESP fit and GAFF parameter generation
  generate_amber_params

  # converting AMBER files to GROMACS
  python $PATH_TO_ACEPYPE -p prmtop -x inpcrd -r > log_acepype.log 2>&1
  gmx editconf -f AA_GMX.gro -o AA_GMX.pdb > log_editconf.log 2>&1

  if [ ! -f AA_GMX.pdb ];then
    echo "automatic ligand parametrization failed in $(pwd), continue?"
    select yn in "Yes" "No"; do
        case $yn in
            Yes ) break;;
            No ) exit;;
        esac
    done
  fi
  echo 'automatic generation of GAFF parameters might have worked...'
  # fix a bug with GROMACS not allowing atom names O1, O2 in residues/ligands
  sed -i 's/O1 /O91/g' $BASEDIR/$FFOUTPATH/$MOLNAME/ligand/AA_GMX.pdb
  sed -i 's/O2 /O92/g' $BASEDIR/$FFOUTPATH/$MOLNAME/ligand/AA_GMX.pdb
  sed -i 's/ O1 / O91/g' $BASEDIR/$FFOUTPATH/$MOLNAME/ligand/AA_GMX.top
  sed -i 's/ O2 / O92/g' $BASEDIR/$FFOUTPATH/$MOLNAME/ligand/AA_GMX.top

  # convert ligand parameters to GROMACS pdb2gmx / grompp format
  mkdir -p $BASEDIR/$FFOUTPATH/$MOLNAME/top2itp
  cd $BASEDIR/$FFOUTPATH/$MOLNAME/top2itp
  cp $BASEDIR/tools/top2itp/* .
  cp $BASEDIR/$FFOUTPATH/$MOLNAME/ligand/AA_GMX.top $BASEDIR/$FFOUTPATH/$MOLNAME/top2itp/topol.top
  python RUNME.py > log_RUNME1.log 2>&1
  python RUNME2.py > log_RUNME2.log 2>&1

  mkdir -p $BASEDIR/$FFOUTPATH/$MOLNAME/gmxVsites
  cd $BASEDIR/$FFOUTPATH/$MOLNAME/gmxVsites
  mkdir -p amber99sb-ildn.ff
  cd amber99sb-ildn.ff

  modify_amber99sb_ff_for_ligand

  cd $BASEDIR/$FFOUTPATH/$MOLNAME/gmxVsites
  cp $BASEDIR/$FFOUTPATH/$MOLNAME/ligand/AA_GMX.pdb mol.pdb

  echo -e '1\n1\n' | gmx pdb2gmx -f mol.pdb -vsite hydrogens -o vs.pdb > log_pdb2gmx.log 2>&1
  cp $BASEDIR/tools/md.mdp .
  gmx editconf -f vs.pdb -o vs_box.pdb -d 1.3 > log_editconf.log 2>&1 || echo 'editconf failed'
  gmx grompp -f md.mdp -c vs_box.pdb -p topol.top -pp processed.top -o top.tpr > log_grompp.log 2>&1 || echo 'grompp failed'; a=$(gmx grompp -f md.mdp -c vs_box.pdb -p topol.top -pp processed.top -o top.tpr 2>&1 | grep Unknown| awk '{print $3}'); echo "attempting to autoremove unkown atom type $a";cat amber99sb-ildn.ff/aminoacids.vsd | egrep -v $a > amber99sb-ildn.ff/aa.vsd; mv amber99sb-ildn.ff/aa.vsd amber99sb-ildn.ff/aminoacids.vsd; cat amber99sb-ildn.ff/ffbonded.itp | egrep -v $a > amber99sb-ildn.ff/aa.bon; mv amber99sb-ildn.ff/aa.bon amber99sb-ildn.ff/ffbonded.itp
  gmx grompp -f md.mdp -c vs_box.pdb -p topol.top -pp processed.top -o top.tpr > log_grompp_rerun.log 2>&1
  if [ ! -f processed.top ];then
    echo "pdb2gmx/grommp generation of V-SITEs failed in $(pwd), continue?"
    select yn in "Yes" "No"; do
        case $yn in
            Yes ) break;;
            No ) exit;;
        esac
    done
  fi
  echo "GAFF V-SITEs generation in $(pwd) seemed to work. ;-)\n\n\t\t FIN \n\n"
  cd $BASEDIR
done
