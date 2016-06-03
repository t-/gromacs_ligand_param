import numpy as np

class merge_dihedrals:
    def __init__(self,filepath,filelist,dihetype='9',fmerged='merged_dihedrals.itp'):
        clist=[]
        for fi in filelist:
            flist1 = self.read_dihedrals(filepath+fi,t=dihetype)
            clist=clist+flist1
            #print 'processing',fi
        print 'fixing type',dihetype,'dihedrals'
        clist=self.fix_permuted_entries(clist)
        clist=self.fix_double_entries(clist)
        clist=self.get_similars(clist,filepath,fname='dihedral_errors.dat')
        #self.print_dihedrals(clist)
        clist.sort()
        self.print2file_dihedrals(clist,filepath,fmerged)

    def print_dihedrals(self,clist):
        for lin in clist:
            clin=lin[0]+' '+lin[1]
            top=clin.split(' ')
            out='%4s%4s%4s%4s %2s%8s%14s%4s' % (top[0],top[1],top[2],top[3],top[4],top[5],top[6],top[7])
            print out

    def print2file_dihedrals(self,clist,fpath,fname):
        f=open(fpath+fname,'w')
        print >>f, '[ dihedraltypes ]'
        print >>f, '; i   j   k   l func  phase     kd      pn'
        for lin in clist:
            clin=lin[0]+' '+lin[1]
            top=clin.split(' ')
            out='%4s%4s%4s%4s %2s%8s%14s%4s' % (top[0],top[1],top[2],top[3],top[4],top[5],top[6],top[7])
            print >>f,out

    def get_similars(self,clist,filepath,fname='dihedral_errors.dat'):
        print 'fixing similar dihedrals - output written to',filepath+fname
        #=======================================================================
        # fixes:
        #  character identical entries with different force constants
        #    cag cag cag cag  9   180.0      15.16700   2
        #    cag cag cag cag  9   180.0      16.73600   2
        # Will always use the larger one of the two by default
        #=======================================================================
        sim_clist={}
        for lin in clist:
            clin=lin[0]+' '+lin[1]
            top=clin.split(' ')
            sim_clist[top[3] + ' ' + top[2] + ' ' + top[1] + ' ' + top[0]+' '+top[4] + ' ' + top[5] + ' ' + top[7]]=[top[6],lin]
        f=open(filepath+fname,'aw')
        print >> f, 'fixed dihedrals'
        for i in xrange(len(clist)):
            lin=clist[i]
            clin=lin[0]+' '+lin[1]
            top=clin.split(' ')
            cur = top[3] + ' ' + top[2] + ' ' + top[1] + ' ' + top[0]+' '+top[4] + ' ' + top[5] + ' ' + top[7]
            if top[6] != sim_clist[cur][0]:
                #This will allways use the larger force constant from the set
                if float(top[6]) > float(sim_clist[cur][0]):
                    print >> f, 'new',top[6],'old',sim_clist[cur][0],sim_clist[cur][1]
                    sim_clist[top[3] + ' ' + top[2] + ' ' + top[1] + ' ' + top[0]+' '+top[4] + ' ' + top[5] + ' ' + top[7]] = [top[6],[top[0] + ' ' + top[1] + ' ' + top[2] + ' ' + top[3],top[4] + ' ' + top[5] + ' ' + top[6] + ' ' + top[7]]]
                if float(top[6]) < float(sim_clist[cur][0]):
                    print >> f, 'new',sim_clist[cur][0],'old',top[6],sim_clist[cur][1]
        new_clist=[]
        f.close()
        for i in sim_clist.keys(): 
            new_clist.append(sim_clist[i][1]) 
        return clist

    def fix_permuted_entries(self,clist):
        print 'fixing permuted dihedrals'
        #=======================================================================
        # fixes:
        #  character identical permuted entries like
        #    nhg c2g ceg hag 9 180.0 27.82360 2
        #    hag ceg c2g nhg 9 180.0 27.82360 2
        #=======================================================================
        perm_clist=[]
        for lin in clist:
            clin=lin[0]+' '+lin[1]
            top=clin.split(' ')
            order=[top[0]+' '+top[1],top[3]+' '+top[2]]
            order_ref=[top[0]+' '+top[1],top[3]+' '+top[2]]
            order_ref.sort()
            if order!=order_ref:
                perm_clist.append([top[3] + ' ' + top[2] + ' ' + top[1] + ' ' + top[0],top[4] + ' ' + top[5] + ' ' + top[6] + ' ' + top[7]])
            else: 
                perm_clist.append([top[0] + ' ' + top[1] + ' ' + top[2] + ' ' + top[3],top[4] + ' ' + top[5] + ' ' + top[6] + ' ' + top[7]])
        return perm_clist

    def fix_double_entries(self,clist):
        print 'fixing double dihedrals'
        #=======================================================================
        # fixes:
        #  character identical entries like
        #    nhg c2g ceg hag 9 180.0 27.82360 2
        #    nhg c2g ceg hag 9 180.0 27.82360 2
        #=======================================================================
        keys = {}
        for e in clist:
            ie=e[0]+' '+e[1]
            keys[ie] = 1
        lins=keys.keys()
        lins.sort()
        #splits list up again and converts it back into input format: ['cag cfg ceg hg','9 180.0 27.82360 2']
        linreturn=[]
        for lin in lins:
            top=lin.split(' ')
            linreturn.append([top[0] + ' ' + top[1] + ' ' + top[2] + ' ' + top[3],top[4] + ' ' + top[5] + ' ' + top[6] + ' ' + top[7]])
        return linreturn 

    def read_dihedrals(self, filename, t='9'):
        if t=='9':
            blockn=3
        if t=='4':
            blockn=2
        block = []
        blocklist = []
        #Read Topology and separate it into blocks [ atoms ], [ bonds ], etc.
        for i in open(filename, 'r'):
            if len(i.strip('\n')) == 0:    # blank line indicates the end of a block [ atoms ], [ bonds ], etc.
                if len(block) > 0: blocklist.append(block);
                block = []
            elif len(i.strip('\n')) > 0:    # read block
                block.append(i.strip('\n'))
        blocklist.append(block);
        dihedralslist = []
        for dihedral in blocklist[blockn]:
            if dihedral[0] != '[' and dihedral[0] != ';':
                top = dihedral.split(' ')
                for i in range(top.count('')): top.remove(''); #remove blanks from array
                dihedralslist.append([top[0] + ' ' + top[1] + ' ' + top[2] + ' ' + top[3],top[4] + ' ' + top[5] + ' ' + top[6] + ' ' + top[7]])
        return dihedralslist

class merge_bonds:
    def __init__(self,filepath,filelist,fmerged='merged_bonds.itp'):
        clist=[]
        for fi in filelist:
            flist1 = self.read_bonds(filepath+fi)
            clist=clist+flist1
            #print 'processing',fi
        clist=self.fix_permuted_entries(clist)
        clist=self.fix_double_entries(clist)
        clist=self.get_similars(clist,filepath,fname='bond_errors.dat')
        #self.print_bonds(clist)
        self.print2file_bonds(clist,filepath,fmerged)

    def print_bonds(self,clist):
        for lin in clist:
            clin=lin[0]+' '+lin[1]
            top=clin.split(' ')
            out='%4s%4s %2s%8s%14s' % (top[0],top[1],top[2],top[3],top[4])
            print out

    def print2file_bonds(self,clist,fpath,fname):
        f=open(fpath+fname,'w')
        print >>f, '[ bondtypes ]'
        print >>f, '; i    j  func       b0          kb'
        for lin in clist:
            clin=lin[0]+' '+lin[1]
            top=clin.split(' ')
            out='%4s%4s %2s%8s%14s' % (top[0],top[1],top[2],top[3],top[4])
            print >>f,out

    def get_similars(self,clist,filepath,fname='bond_errors.dat'):
        print 'fixing similar bonds - output written to',filepath+fname
        #=======================================================================
        # fixes:
        #  character identical entries with different force constants
        #     cag cag  1  0.1387      400330.0
        #     cag cag  1  0.1429      350030.0
        #=======================================================================
        sim_clist={}
        for lin in clist:
            clin=lin[0]+' '+lin[1]
            top=clin.split(' ')
            sim_clist[top[0] + ' ' + top[1]]=[top[2] + ' ' + top[3] + ' ' + top[4],[lin]]

        for lin in clist:
            clin=lin[0]+' '+lin[1]
            top=clin.split(' ')
            cur = top[0] + ' ' + top[1]
            if top[2] + ' ' + top[3] + ' ' + top[4] != sim_clist[cur][0]:
                sim_clist[cur][1].append([top[0] + ' ' + top[1],top[2] + ' ' + top[3] + ' ' + top[4]])
        f=open(filepath+fname,'w')
        for lin in sim_clist.keys():
            dmean=[]
            kmean=[]
            if len(sim_clist[lin][1])>1:
                for element in sim_clist[lin][1]: 
                    dmean.append(float(element[1].split(' ')[1]))
                    kmean.append(float(element[1].split(' ')[2]))
                print >>f,'\nBOND TYPE  ',sim_clist[lin][1][0][0]
                print >>f,' distances ',np.array(dmean)
                print >>f,'  mean',np.array(dmean).mean(),'+\-',np.array(dmean).std()
                print >>f,' forceconstants',np.array(kmean)
                print >>f,'  mean',np.array(kmean).mean(),'+\-',np.array(kmean).std()
                #replacing old bond with new averaged bond parameters
                sim_clist[lin][0] = '1 '+str(np.round(np.array(dmean).mean(),4))+' '+str(np.round(np.array(kmean).mean(),0))
        f.close()
        #creating new clist with averaged bond parameters
        new_clist=[]
        for i in sim_clist.keys():
            new_clist.append([i,sim_clist[i][0]])
        new_clist.sort()
        return new_clist

    def fix_permuted_entries(self,clist):
        print 'fixing permuted bonds'
        #=======================================================================
        # fixes:
        #  character identical permuted entries like
        #    cag osg        1     0.1373   311620.0
        #    osg cag        1     0.1373   311620.0
        #=======================================================================
        perm_clist=[]
        for lin in clist:
            clin=lin[0]+' '+lin[1]
            top=clin.split(' ')
            order=[top[0],top[1]]
            order_ref=[top[1],top[0]]
            order_ref.sort()
            if order!=order_ref:
                perm_clist.append([top[1] + ' ' + top[0],top[2] + ' ' + top[3] + ' ' + top[4]])
            else: 
                perm_clist.append([top[0] + ' ' + top[1],top[2] + ' ' + top[3] + ' ' + top[4]])
        return perm_clist

    def fix_double_entries(self,clist):
        print 'fixing double bonds'
        #=======================================================================
        # fixes:
        #  character identical entries like
        #     cag cag        1     0.1429   350030.0 
        #     cag cag        1     0.1429   350030.0 
        #=======================================================================
        keys = {}
        for e in clist:
            ie=e[0]+' '+e[1]
            keys[ie] = 1
        lins=keys.keys()
        lins.sort()
        #splits list up again and converts it back into input format: ['cag cfg ceg hg','9 180.0 27.82360 2']
        linreturn=[]
        for lin in lins:
            top=lin.split(' ')
            linreturn.append([top[0] + ' ' + top[1],top[2] + ' ' + top[3] + ' ' + top[4]])
        return linreturn 

    def read_bonds(self, filename):
        block = []
        blocklist = []
        #Read Topology and separate it into blocks [ atoms ], [ bonds ], etc.
        for i in open(filename, 'r'):
            if len(i.strip('\n')) == 0:    # blank line indicates the end of a block [ atoms ], [ bonds ], etc.
                if len(block) > 0: blocklist.append(block);
                block = []
            elif len(i.strip('\n')) > 0:    # read block
                block.append(i.strip('\n'))
        blocklist.append(block);
        bondslist = []
        for bond in blocklist[0]:
            if bond[0] != '[' and bond[0] != ';':
                top = bond.split(' ')
                for i in range(top.count('')): top.remove(''); #remove blanks from array
                bondslist.append([top[0] + ' ' + top[1],top[2] + ' ' + top[3] + ' ' + top[4]])
        return bondslist

class merge_angles:
    def __init__(self,filepath,filelist,fmerged='merged_angles.itp'):
        clist=[]
        for fi in filelist:
            flist1 = self.read_angles(filepath+fi)
            clist=clist+flist1
            #print 'processing',fi
        clist=self.fix_permuted_entries(clist)
        clist=self.fix_double_entries(clist)
        clist.sort()
        clist=self.get_similars(clist,filepath,fname='angle_errors.dat')
        #self.print_angles(clist)
        self.print2file_angles(clist,filepath,fmerged)

    def fix_permuted_entries(self,clist):
        print 'fixing permuted angles'
        #=======================================================================
        # fixes:
        #  character identical permuted entries like
        #    ssg c3g h1g          1   109.340    449.030 ; TTT
        #    h1g c3g ssg          1   109.340    449.030 ; TTT
        #=======================================================================
        perm_clist=[]
        for lin in clist:
            clin=lin[0]+' '+lin[1]
            top=clin.split(' ')
            order=[top[0],top[2]]
            order_ref=[top[2],top[0]]
            order_ref.sort()
            if order!=order_ref:
                perm_clist.append([top[2] + ' ' + top[1] + ' ' + top[0], top[3] + ' ' + top[4] + ' ' + top[5]])
            else: 
                perm_clist.append([top[0] + ' ' + top[1] + ' ' + top[2], top[3] + ' ' + top[4] + ' ' + top[5]])
        return perm_clist

    def fix_double_entries(self,clist):
        print 'fixing double angles'
        #=======================================================================
        # fixes:
        #  character identical entries like
        #    ssg c3g h1g          1   109.340    449.030 ; TTT
        #    ssg c3g h1g          1   109.340    449.030 ; TTT
        #=======================================================================
        keys = {}
        for e in clist:
            ie=e[0]+' '+e[1]
            keys[ie] = 1
        lins=keys.keys()
        lins.sort()
        #splits list up again and converts it back into input format: ['cag cfg ceg','9 180.0 27.82360']
        linreturn=[]
        for lin in lins:
            top=lin.split(' ')
            linreturn.append([top[0] + ' ' + top[1] + ' ' + top[2], top[3] + ' ' + top[4] + ' ' + top[5]])
        return linreturn 

    def get_similars(self,clist,filepath,fname='angle_errors.dat'):
        print 'fixing similar angles - output written to',filepath+fname
        #=======================================================================
        # fixes:
        #  character identical entries with different force constants
        #     ssg c3g h1g          1   109.340    449.030 ; TTT
        #     ssg c3g h1g          1   29.340     142.030 ; TTT
        #=======================================================================
        sim_clist={}
        for lin in clist:
            clin=lin[0]+' '+lin[1]
            top=clin.split(' ')
            sim_clist[top[0] + ' ' + top[1] + ' ' + top[2]]=[top[3] + ' ' + top[4] + ' ' + top[5],[lin]]

        for lin in clist:
            clin=lin[0]+' '+lin[1]
            top=clin.split(' ')
            cur = top[0] + ' ' + top[1] + ' ' + top[2]
            if top[3] + ' ' + top[4] + ' ' + top[5] != sim_clist[cur][0]:
                sim_clist[cur][1].append([top[0] + ' ' + top[1] + ' ' + top[2], top[3] + ' ' + top[4] + ' ' + top[5]])
        f=open(filepath+fname,'w')
        for lin in sim_clist.keys():
            dmean=[]
            kmean=[]
            if len(sim_clist[lin][1])>1:
                for element in sim_clist[lin][1]: 
                    dmean.append(float(element[1].split(' ')[1]))
                    kmean.append(float(element[1].split(' ')[2]))
                print >>f,'\nAngle TYPE  ',sim_clist[lin][1][0][0]
                print >>f,' distances ',np.array(dmean)
                print >>f,'  mean',np.array(dmean).mean(),'+\-',np.array(dmean).std()
                print >>f,' forceconstants',np.array(kmean)
                print >>f,'  mean',np.array(kmean).mean(),'+\-',np.array(kmean).std()
                #replacing old bond with new averaged bond parameters
                sim_clist[lin][0] = '1 '+str(np.round(np.array(dmean).mean(),4))+' '+str(np.round(np.array(kmean).mean(),0))
        f.close()
        #creating new clist with averaged bond parameters
        new_clist=[]
        for i in sim_clist.keys():
            new_clist.append([i,sim_clist[i][0]])
        new_clist.sort()
        return new_clist

    def read_angles(self, filename):
        block = []
        blocklist = []
        #Read Topology and separate it into blocks [ atoms ], [ bonds ], etc.
        for i in open(filename, 'r'):
            if len(i.strip('\n')) == 0:    # blank line indicates the end of a block [ atoms ], [ bonds ], etc.
                if len(block) > 0: blocklist.append(block);
                block = []
            elif len(i.strip('\n')) > 0:    # read block
                block.append(i.strip('\n'))
        blocklist.append(block);
        angleslist = []
        for angle in blocklist[1]:
            if angle[0] != '[' and angle[0] != ';':
                top = angle.split(' ')
                for i in range(top.count('')): top.remove(''); #remove blanks from array
                angleslist.append([top[0] + ' ' + top[1] + ' ' + top[2] , top[3] + ' ' + top[4] + ' ' +top[5]])
        return angleslist
    
    def print_angles(self,clist):
        for lin in clist:
            clin=lin[0]+' '+lin[1]
            top=clin.split(' ')
            out='%4s%4s%4s %8s%14s%14s' % (top[0],top[1],top[2],top[3],top[4],top[5])
            print out

    def print2file_angles(self,clist,fpath,fname):
        f=open(fpath+fname,'w')
        print >>f, '[ angletypes ]'
        print >>f, '; i    j    k   func       phi0          k'
        for lin in clist:
            clin=lin[0]+' '+lin[1]
            top=clin.split(' ')
            out='%4s%4s%4s %8s%14s%14s' % (top[0],top[1],top[2],top[3],top[4],top[5])
            print >>f,out

class merge_atomtypes:
    def __init__(self,filepath,filelist,fmerged='merged_atomtypes.atp'):
        clist=[]
        for fi in filelist:
            flist1 = self.read_atomtypes(filepath+fi)
            clist=clist+flist1
        clist=self.fix_double_entries(clist)
        clist.sort()
        self.print2file_angles(clist,filepath,fmerged)

    def read_atomtypes(self,filename):
        clist=[]
        for i in open(filename):
            clist.append( i.strip('\n') )
        return clist
    
    def fix_double_entries(self,clist):
        print 'fixing double atomtypes'
        #=======================================================================
        # fixes:
        #  character identical entries like
        #   n2g               14.01000    ; TTT
        #   n2g               14.01000    ; TTT
        #=======================================================================
        keys = {}
        for e in clist:
            keys[e] = 1
        lins=keys.keys()
        lins.sort()
        #splits list up again and converts it back into input format: ['cag cfg ceg','9 180.0 27.82360']
        return lins 
    
    def print2file_angles(self,clist,fpath,fname):
        f=open(fpath+fname,'w')
        for lin in clist:
            print >>f,lin

class merge_nonbonded:
    def __init__(self,filepath,filelist,fmerged='merged_nbonds.itp'):
        clist=[]
        for fi in filelist:
            flist1 = self.read_atomtypes(filepath+fi)
            clist=clist+flist1
            #print 'processing',fi
        clist=self.fix_double_entries(clist)
        clist.sort()
        self.print2file_angles(clist,filepath,fmerged)

    def read_atomtypes(self,filename):
        clist=[]
        for i in open(filename):
            if i.find('[')<0:
                clist.append( i.strip('\n') )
        return clist

    def fix_double_entries(self,clist):
        print 'fixing double nonbonded parameters'
        #=======================================================================
        # fixes:
        #  character identical entries like
        #   ohg          8       16.0    0.0000  A  3.066470e-01 8.803140e-01
        #   ohg          8       16.0    0.0000  A  3.066470e-01 8.803140e-01
        #=======================================================================
        keys = {}
        for e in clist:
            keys[e] = 1
        lins=keys.keys()
        lins.sort()
        return lins 

    def print2file_angles(self,clist,fpath,fname):
        f=open(fpath+fname,'w')
        print >>f, '[ atomtypes ]'
        for lin in clist:
            print >>f,lin

def main():
    fpath='./'
    print 'working in directory',fpath
    f=open(fpath+'dihedral_errors.dat','w')
    print >>f,''
    f.close()
    merge_dihedrals('./',['res_ffbonded.itp'],dihetype='9',fmerged='merged_dihedrals.itp')
    print ''
    merge_dihedrals('./',['res_ffbonded.itp'],dihetype='4',fmerged='merged_impropers.itp')
    print ''
    merge_bonds('./',['res_ffbonded.itp'],fmerged='merged_bonds.itp')
    print ''
    merge_angles('./',['res_ffbonded.itp'],fmerged='merged_angles.itp')
    print ''
    merge_atomtypes('./',['res_atomtypes.atp'],fmerged='merged_atomtypes.atp')
    print ''
    merge_nonbonded('./',['res_ffnonbonded.itp'],fmerged='merged_nbonds.itp')

if __name__ == '__main__':
    main()
