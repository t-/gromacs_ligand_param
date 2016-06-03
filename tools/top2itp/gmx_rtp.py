class gmx_rtp():
    def __init__(self):
        self.resname     = 'MOL'
        self.atoms       = []
        self.bonds       = []
        self.angles      = []
        self.dihedrals   = []
        self.impropers   = []
    def rtp_gen(self,top):
        def atype(number,top):
            ###gets an atom number and returns its atom name from the [ atoms ] section in the .top file (top object ->top.atoms_list)
            ai = (int)(number)
            for i in top.atoms_list:
                if i.nr == ai:
                    return i.atom       # 'atom' = atomname as found in the .top file header 
            return 'XX'                 #default value if nothing is found
        def isa(n,top):#check if n is amber atom type=1 or gaff=0
            ai = (int)(n)
            for i in top.atoms_list:
                if i.nr == ai:
                    if i.type.find('amber')>=0:
                        return 1
                    else:
                        return 0                    
        for at in top.atoms_list:       #define all the atoms in ffamber99sb.rtp format
            tmp = self.atom()
            tmp.name   = at.atom
            tmp.type   = at.type
            tmp.charge = at.charge
            self.atoms.append(tmp)
            for i,atom in enumerate(self.atoms):
                atom.number=i+1
        for bn in top.bonds_list:       #define all the bonds in ffamber99sb.rtp format 
	    tmp = self.bond()
            tmp.i = atype(bn.ai,top)
            tmp.j = atype(bn.aj,top)
            tmp.r = bn.r
            tmp.k = bn.k
            if not (isa(bn.ai,top)==1 and isa(bn.aj,top)==1): #dont add Amber/Amber bonds, only Amber/gaff and gaff/gaff
                self.bonds.append(tmp)
            else:
                tmp.r=''
                tmp.k=''
                self.bonds.append(tmp)
        #self.bonds.append(self.bond('-C','N'))              #this is needed to connect our residue to the next one :-)
        '''
        for an in top.angles_list:
            tmp = self.angle()
            tmp.ai = atype(an.ai,top)
            tmp.aj = atype(an.aj,top)
            tmp.ak = atype(an.ak,top)
            tmp.theta = an.theta 
            tmp.cth   = an.cth
            if not (isa(an.ai,top)==1 and isa(an.aj,top)==1 and isa(an.ak,top)==1): #dont add Amber/Amber bonds, only Amber/gaff and gaff/gaff
                self.angles.append(tmp)
         
        for di in top.dihedrals_list:
            tmp = self.dihedral()
            tmp.i = atype(di.i,top)
            tmp.j = atype(di.j,top)
            tmp.k = atype(di.k,top)
            tmp.l = atype(di.l,top)
            tmp.C0= di.phase
            tmp.C1= di.kd
            tmp.C2= di.pn
            if not (isa(di.i,top)==1 and isa(di.j,top)==1 and isa(di.k,top)==1 and isa(di.l,top)==1): #dont add Amber/Amber bonds, only Amber/gaff and gaff/gaff
                self.dihedrals.append(tmp)
	'''
        for di in top.impropers_list:
            tmp = self.improper()
            tmp.i = atype(di.i,top)
            tmp.j = atype(di.j,top)
            tmp.k = atype(di.k,top)
            tmp.l = atype(di.l,top)
            tmp.C0= di.phase
            tmp.C1= di.kd
            tmp.C2= di.pn
            if not (isa(di.i,top)==1 and isa(di.j,top)==1 and isa(di.k,top)==1 and isa(di.l,top)==1): #dont add Amber/Amber bonds, only Amber/gaff and gaff/gaff
                self.impropers.append(tmp)


        #print 'generating connecting backbone dihedrals!'   #same for dihedrals 
        #self.dihedrals.append(self.dihedral( 'CA'  ,   'C'  ,  '+N'  ,  '+H'  ,  'backbone_prop_1'))
        #self.dihedrals.append(self.dihedral(  'O'  ,   'C'  ,  '+N'  ,  '+H'  ,  'backbone_prop_2'))
        #self.dihedrals.append(self.dihedral( '-C'  ,   'N'  ,  'CA'  ,  'CB'  ,  'backbone_prop_3'))
        #self.dihedrals.append(self.dihedral( '-C'  ,   'N'  ,  'CA'  ,   'C'  ,  'backbone_prop_4'))
        #self.dihedrals.append(self.dihedral( 'CA'  ,   'C'  ,  '+N'  , '+CA'  ,  'backbone_prop_1'))
        #self.dihedrals.append(self.dihedral(  'O'  ,   'C'  ,  '+N'  , '+CA'  ,  'backbone_prop_1'))
        #print 'generating connecting backbone impropers!'   #same for impropers
        #self.impropers.append(self.improper( '-C','CA','N','H' ))
        #self.impropers.append(self.improper( 'CA','+N','C','O' ))
        
    def rtp_print(self,filename,top):
        print 'writing',filename
        #use the cprint() method of each object to produce the proper rtp file format (see below)! 
        f = open(filename,'w')
        print >>f, '[ '+self.resname+' ]'
        print >>f, ' [ atoms ]'
        for i in self.atoms:
            print >> f, i.cprint()
        print >>f, ' [ bonds ]'
        for x in self.bonds:
            #print x.cprint()
            print >> f,x.c45print()
        '''
        print >>f,' [ angles ] '
        for i in self.angles:
            #print i.cprint()
            print >>f, i.c45print()
        print >>f, ' [ dihedrals ]'
        for i in self.dihedrals:
            print >> f,i.c45print()
	'''
        print >>f, ' [ impropers ]'
        for i in self.impropers:
            print >> f,i.c45print()
        f.close()

    class atom():
        def __init__(self,name='',type='',charge=0.0,number=0):
            self.name = str(name)
            self.type = str(type)
            try:
                self.charge=(float)(charge)
            except Exception:
                self.charge=0.0
            try:
                self.number=number
            except Exception:
                self.number=0
        def r(self,string,length):
            string = str(string)
            while len(string)<length: string= ' ' + string;
            return string
        def l(self,string,length):
            string = str(string)
            while len(string)<length: string= string+' ';
            return string
        def cprint(self):
            return self.r(self.name,6)+'    '+self.l(self.type,10)+self.r(self.charge,10)+self.r(self.number,6)

    class bond():
        def __init__(self,i='',j='',r='',k=''):
            self.i = str(i)
            self.j = str(j)
            self.r = str(r)
            self.k = str(k)
        def s(self,string,length):
            string = str(string)
            while len(string)<length: string= ' ' + string;
            return string
        def cprint(self):
            return self.s(self.i,6) + self.s(self.j,6)+ self.s(1,6) + self.s(self.r,12) + self.s(self.k,12)
        def c45print(self):
            return self.s(self.i,6) + self.s(self.j,6)

    class angle():
        def __init__(self, ai='',aj='',ak='',theta='',cth=''):
            self.ai      = ai.strip(' ')
            self.aj      = aj.strip(' ')
            self.ak      = ak.strip(' ')
            try:
                self.theta   = theta.strip(' ')
                self.cth     = cth.strip(' ')
            except Exception:
                self.theta   = ''
                self.cth     = ''
        def r(self,string,length):
            string = str(string)
            while len(string)<length: string= ' ' + string;
            return string
        def cprint(self):
            return self.r(self.ai,5) + self.r(self.aj,6) + self.r(self.ak,6)+ self.r(1,6) + self.r(self.theta,12) + self.r(self.cth,12)
        def c45print(self):
	    return self.r(self.ai,5) + self.r(self.aj,6) + self.r(self.ak,6)
    class dihedral():
        def __init__(self,i='',j='',k='',l='',C0='',C1='',C2=''):
            self.i = str(i)
            self.j = str(j)
            self.k = str(k)
            self.l = str(l)
            self.C0= str(C0)
            self.C1= str(C1)
            self.C2= str(C2)
                    

        def r(self,string,length):
            string = str(string)
            while len(string)<length: string= ' ' + string;
            return string
        def cprint(self):
            return self.r(self.i,6)+self.r(self.j,6)+self.r(self.k,6)+self.r(self.l,6)+self.r(9,6)+self.r(self.C0,16)+self.r(self.C1,10)+self.r(self.C2,10)
        def c45print(self):
            return self.r(self.i,6)+self.r(self.j,6)+self.r(self.k,6)+self.r(self.l,6)
    	
    class improper():
        def __init__(self,i='',j='',k='',l='',type='',C0='',C1='',C2=''):
            self.i = str(i)
            self.j = str(j)
            self.k = str(k)
            self.l = str(l)
            self.C0= str(C0)
            self.C1= str(C1)
            self.C2= str(C2)
            self.type=str(type)
        def r(self,string,length):
            string = str(string)
            while len(string)<length: string= ' ' + string;
            return string
        def cprint(self):
            return self.r(self.i,6)+self.r(self.j,6)+self.r(self.k,6)+self.r(self.l,6)+self.r(9,6)+self.r(self.C0,16)+self.r(self.C1,10)+self.r(self.C2,10)
        def c45print(self):
            return self.r(self.i,6)+self.r(self.j,6)+self.r(self.k,6)+self.r(self.l,6)
