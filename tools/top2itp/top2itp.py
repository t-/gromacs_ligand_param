from gmx_top import *
from gmx_bon import *
from gmx_nb import *
from gmx_rtp import *
from gmx_atp import *
import os

class top2itp():
    def __init__(self,topfilename,nbfile='res_ffnonbonded.itp',bonfile='res_ffbonded.itp',atpfile='res_atomtypes.atp',rtpfile='res_aminoacids.rtp'):
        self.top =topology()                #GMX top object
        self.top.read_top(topfilename)      #read in the gaff residue.top file from amb2gmx.pl
        #Fixing CW ring atoms
        cw1=os.environ.get('cw1')
        cw2=os.environ.get('cw2')
        if (cw1!=None) and (cw2!=None):
            print 'Fixing CW ring atoms: cw1',cw1,' cw2',cw2
            a=self.top.atomtype('cw','cw','0.00000','0.00000','A','3.39967e-01','3.59824e-01')
            print a.cprint()
            self.top.atomtypes_list.append(a)
            for i in self.top.atoms_list:
                if (cw1==i.atom) or (cw2==i.atom):
                     print i.atom,cw1,cw2
                     i.type='cw'
        #/Fixing CW ring atoms

        #Fix CA Ring
        calist=os.environ.get('calist')
        if (calist!=None):
            casplit=calist.split(',')
            for i in casplit:
                print 'Fixing CA ring atoms:', i
            for i in self.top.atoms_list:
                if i.atom in casplit:
                     print i.atom,casplit
                     i.type='ca'
        #/Fix CA Ring

        #Fix C2 - C Atoms
        c2=os.environ.get('c2')
        if (c2!=None):
            print 'Fixing C2 Atom ',c2
            for i in self.top.atoms_list:
                if i.atom == c2:
                     print i.atom,casplit
                     i.type='c'
        #/Fix C2 - C Atoms
        
        #Fixing C O CA HA N H
        link_type=os.environ.get('LINK_TYPE')
        link_n=os.environ.get('LINK_N')
        link_h=os.environ.get('LINK_H')
        link_ca=os.environ.get('LINK_CA')
        link_ha=os.environ.get('LINK_HA')
        link_c=os.environ.get('LINK_C')
        link_o=os.environ.get('LINK_O')
        print link_type
        print link_n
        print link_h
        print link_ca
        print link_ha
        print link_c
        print link_o
                                                
        for i in self.top.atoms_list:        
            if (link_n!=None) and (link_n==i.atom):
                i.atom='N'
                i.type='N'
            if (link_h!=None) and (link_h==i.atom):
                i.atom='H'
                i.type='H'
            if (link_ca!=None) and (link_ca==i.atom):
                i.atom='CA'
                i.type='CT'
            if (link_ha!=None) and (link_ha==i.atom):
                i.atom='HA'
                i.type='H1'
                if (link_type=='Nter'):
                    i.type='HP'
            if (link_c!=None) and (link_c==i.atom):
                i.atom='C'
                i.type='C'                
            if (link_o!=None) and (link_o==i.atom):
                i.atom='O'
                i.type='O'                
        #/Fixing C O CA HA N H        

        #Fix O1, O2 issues
        for i in self.top.atoms_list:
            if (i.atom=='O1'):
                i.atom='O91'
            if (i.atom=='O2'):
                i.atom='O92'
        #Fix O1, O2 issues

        for atm in self.top.atoms_list:
            if atm.type not in ['C','H','N','HC','H1','O','O2','N3','CT','HP']:
                atm.type=atm.type[0]+atm.type[1:]+'g'
        for atm in self.top.atomtypes_list:
            if atm.name not in ['C','H','N','HC','H1','O','O2','N3','CT','HP']:
                atm.name=atm.name[0]+atm.name[1:]+'g'
        self.nb  = ffnb()                   #object for non-bonded parameters
        self.top2nb(self.top)               #convert gmx.top to ffamber99sbnb.itp format
        self.writenb(nbfile)                #write the nb.itp file
        self.bon = ffbon()                  #object for bonded parameters
        self.top2bon(self.top)              #convert gmx.top to ffamber99sbbon.itp format

        #/Fixing CA ring atoms
        #calist=os.environ.get('calist')
        #if (calist!=None):
        #    casplit=calist.split(',')
        #    for i in casplit:
        #        print 'Fixing CA ring atoms:', i
        #    for i in self.top.atoms_list:
        #        if i.atom in casplit:
        #             print i.atom,casplit
        #             i.type='cag'
        #/Fixing CA ring atoms

        ##Fixing CW ring atoms
        for ii in self.bon.dihedrals_list:
            if (ii.j=='cwg') and (ii.k=='cwg'):
                ii.phase=0.0
        for ii in self.bon.impropers_list:
            if (ii.j=='cwg') and (ii.k=='cwg'):
                ii.phase=0.0
        ##/Fixing CW ring atoms

        self.writebon(bonfile)              #write the bon.itp file
        self.top2rtp(self.top,rtpfile)      #takes self.top and converts it to ffamber99sb.rtp file format
        self.atp = gmx_atp(self.top)        #creates the atp object 
        self.atp.write_atp(atpfile)         #writes the atp file in ffamber99sb.atp format
        
    def top2nb(self,top):
        #########################################################
        # will fill the non-bonded parameter object (self.nb)   #
        #########################################################
        for at in top.atomtypes_list:
            tmp = self.nb.atomtype(at.name,str(0.0),str(at.mass),str(at.charge),str(at.ptype),str(at.sigma),str(at.epsilon))
            type= at.name[:1]
            if type == 'c':
                tmp.atnum=6
                tmp.mass=12.01
            if type == 'n':
                tmp.atnum=7
                tmp.mass=14.01
            if type == 'h':
                tmp.atnum=1
                tmp.mass=1.008
            if type == 'o':
                tmp.atnum=8
                tmp.mass=16.00
            if type == 's':
                tmp.atnum=16
                tmp.mass=32.06
            if type == 'f':
                tmp.atnum=9
                tmp.mass=19.00
            if type == 'b':
                tmp.atnum=5
                tmp.mass=10.81
            if type == 'p':
                tmp.atnum=15
                tmp.mass=30.973
            if type not in ['c','n','o','h','s','f','b','p']:
                print 'ERROR-- atomtype',type,' not yet supported here !'
            self.nb.atomtypes.append(tmp)
    
    def top2bon(self,top):
        self.bon.top2bon(top)
        
    def top2rtp(self,top,rtpfile):
        rtp = gmx_rtp()
        rtp.rtp_gen(top)
        rtp.rtp_print(rtpfile,top)

    def writenb(self,nbfile):
        self.nb.write(nbfile)
    
    def writebon(self,bonfile):
        self.bon.write(bonfile)

