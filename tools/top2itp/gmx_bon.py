def s(string,length):
    string = str(string)
    while len(string)<length: string=' ' + string;
    #if len(string)>length: string=string[0:length];
    return string

class ffbon:
    def __init__(self):  
        self.top                 = ''
        self.bonds_list          = []
        self.pairs_list          = []
        self.angles_list         = []
        self.dihedrals_list      = []
        self.impropers_list      = []

        self.pairs_meta         = []
        self.bonds_meta         = []
        self.angles_meta        = []
        self.dihedrals_meta     = []
        self.impropers_meta     = []
    def write(self,filename):
        f = open(filename,'w')
        print 'writing',filename
        print >>f, '[ bondtypes ]'
        print >>f,'; i    j  func       b0          kb'
        for i in self.bonds_list:
            print >> f, i.cprint()
        print >>f,'\n[ angletypes ]'
        print >>f,';  i    j    k  func       th0       cth'
        for i in self.angles_list:
            print >> f, i.cprint()
        print >>f,'\n[ dihedraltypes ]'
        print >>f,'; i   j   k   l func  phase     kd      pn'
        for i in self.impropers_list:
            print >> f, i.cprint()
        print >>f,'\n[ dihedraltypes ]'
        print >>f,'; i   j   k   l func  phase     kd      pn'
        for i in self.dihedrals_list:
            print >> f, i.cprint()
        f.close()
    def type(self,number):
        #here we need to assign a amber99sb bondtype to amber99sb atom types :-)
        ai = (int)(number)
        amber = [['amber99_1','BR'],['amber99_2','C '],['amber99_3','CA'],['amber99_4','CB'],['amber99_5','CC'],['amber99_6','CK'],['amber99_7','CM'],['amber99_8','CN'],['amber99_9','CQ'],['amber99_10','CR'],['amber99_11','CT'],['amber99_12','CV'],['amber99_13','CW'],['amber99_14','C*'],['amber99_15','Ca'],['amber99_16','F '],['amber99_17','H '],['amber99_18','HC'],['amber99_19','H1'],['amber99_20','H2'],['amber99_21','H3'],['amber99_22','HA'],['amber99_23','H4'],['amber99_24','H5'],['amber99_25','HO'],['amber99_26','HS'],['amber99_27','HW'],['amber99_28','HP'],['amber99_29','I '],['amber99_30','Cl'],['amber99_31','Na'],['amber99_32','IB'],['amber99_33','Mg'],['amber99_34','N '],['amber99_35','NA'],['amber99_36','NB'],['amber99_37','NC'],['amber99_38','N2'],['amber99_39','N3'],['amber99_40','N*'],['amber99_41','O '],['amber99_42','OW'],['amber99_43','OH'],['amber99_44','OS'],['amber99_45','O2'],['amber99_46','P '],['amber99_47','S '],['amber99_48','SH'],['amber99_49','CU'],['amber99_50','FE'],['amber99_51','K '],['amber99_52','Rb'],['amber99_53','Cs'],['amber99_54','OW'],['amber99_55','HW'],['amber99_56','Li'],['amber99_57','Zn'],['amber99_58','Sr'],['amber99_59','Ba'],['amber99_60','HW'],['amber99_61','OW'],['amber99_62','HW'],['amber99_63','OW'],['amber99_64','HW'],['amber99_65','OW'],['amber99_66','HW'],['amber99_67','OW'],['CT','CT'],['H1','H1'],['C','C'],['O','O'],['N','N'],['N3','N3'],['HC','HC'],]
        for i in self.top.atoms_list:
            if i.nr == ai:  #not an amber99sb atom type
                if i.type.find('amber')<0:
                    return i.type
                else:       #this is an amber99sb atom type
                    for type in amber:
                        #print type[0]
                        if type[0]==i.type:
                            return type[1]
        return 'XX'
        
    def top2bon(self,top):
        def caps(string):
            ref = string.upper()
            if string == ref:
                return 1
            return 0
        
        self.top = top
        #bonds
        tmpbonds = {}
        for i in top.bonds_list:
            tmpbonds[(self.type(i.ai),self.type(i.aj),str(i.funct),str(i.r),str(i.k))]=''
            #tmp = self.bonds(self.type(i.ai),self.type(i.aj),str(i.funct),str(i.r),str(i.k))
        for i in tmpbonds.keys():
            if caps(i[0]) + caps(i[1]) != 2:  # sort out N - C or CA - C parameters --> those are already in amber, yet keep N - c1 or ha - C 
                #print i[0],i[1]
                self.bonds_list.append(self.bonds(i[0],i[1],i[2],i[3],i[4]))
        #angles
        tmpangles = {}
        for i in top.angles_list:
            tmpangles[(self.type(i.ai),self.type(i.aj),self.type(i.ak),str(i.funct),str(i.theta),str(i.cth))]=''
            #tmp = self.angles(self.type(i.ai),self.type(i.aj),self.type(i.ak),str(i.funct),str(i.theta),str(i.cth))
        for i in tmpangles.keys():
            if caps(i[0]) + caps(i[1]) + caps(i[2])!=3:
                #print i[0],i[1],i[2]
                self.angles_list.append(self.angles(i[0],i[1],i[2],i[3],i[4],i[5]))
        #dihedrals
        tmpdihe = {}
        for x in top.dihedrals_list:
            tmpdihe[(self.type(x.i),self.type(x.j),self.type(x.k),self.type(x.l),str(x.funct),str(x.phase),str(x.kd),str(x.pn))]=''
            #tmp = self.dihedrals(self.type(x.i),self.type(x.j),self.type(x.k),self.type(x.l),str(x.funct),str(x.C0),str(x.C1),str(x.C2),str(x.C3),str(x.C4),str(x.C5))
        for i in tmpdihe.keys():
            if caps(i[0]) + caps(i[1]) + caps(i[2]) + caps(i[3]) != 4:
                #print i[0],i[1],i[2],i[3]
                self.dihedrals_list.append(self.dihedrals(i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7]))
        #impropers
        tmpimp = {}
        for x in top.impropers_list:
            tmpimp[(self.type(x.i),self.type(x.j),self.type(x.k),self.type(x.l),str(x.funct),str(x.phase),str(x.kd),str(x.pn))]=''
        for i in tmpimp.keys():
            if caps(i[0]) + caps(i[1]) + caps(i[2]) + caps(i[3]) != 4:
                self.impropers_list.append(self.impropers(i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7]))
    class bonds:
        def __init__(self,ai,aj,funct,b0='',kb=''):
            self.ai      = (ai.strip(' '))
            self.aj      = (aj.strip(' '))
            self.funct   = (int)(funct.strip(' '))
            try:
                self.b0       = (float)(b0.strip(' '))
                self.kb       = (float)(kb.strip(' '))
            except Exception:
                self.b0       = ''
                self.kb       = ''
        def s(self,string,length):
            string = str(string)
            while len(string)<length: string=' ' + string;
            #if len(string)>length: string=string[0:length];
            return string
        def l(self,string,length):
            string = str(string)
            while len(string)<length: string= string + ' ';
            #if len(string)>length: string=string[0:length];
            return string
        def cprint(self):
            return self.s('',2)+self.l(str(self.ai),3) +self.s('',1)+ self.l(str(self.aj),11)+str(self.funct)+self.s(str('%1.4f'%self.b0),11)+self.s(str('%6.1f'%self.kb),11)+' ; TTT'
    
    class angles:
        #;  i    j    k  func       th0       cth
        #HW  OW  HW           1   104.520    836.800 ; TIP3P water
        def __init__(self, ai,aj,ak,func,th0='',cth=''):
            self.ai      = (ai.strip(' '))
            self.aj      = (aj.strip(' '))
            self.ak      = (ak.strip(' '))
            self.funct   = (int)(func.strip(' '))
            try:
                self.th0     = (float)(th0.strip(' '))
                self.cth     = (float)(cth.strip(' '))
            except Exception:
                self.th0   = ''
                self.cth     = ''
        def s(self,string,length):
            string = str(string)
            while len(string)<length: string=' ' + string;
            #if len(string)>length: string=string[0:length];
            return string
        def l(self,string,length):
            string = str(string)
            while len(string)<length: string= string + ' ';
            #if len(string)>length: string=string[0:length];
            return string
        def cprint(self):
            return self.l(str(self.ai),4)+self.l(str(self.aj),4)+self.l(str(self.ak),13)+self.l(str(self.funct),1)+self.s(str('%3.3f'%self.th0),10)+self.s('%3.3f'%self.cth,11)+' ; TTT'
    
    class dihedrals:
        def __init__(self,i,j,k,l,funct,phase='',kd='',pn=''):
            self.i       = (i.strip(' '))
            self.j       = (j.strip(' '))
            self.k       = (k.strip(' '))
            self.l       = (l.strip(' '))
            self.funct   = (int)(funct.strip(' '))
            try:
                self.phase   = (float)(phase.strip(' '))
                self.kd      = (float)(kd.strip(' '))
                self.pn      = (float)(pn.strip(' ').strip('\t;').strip('\n'))
            except Exception:
                self.phase   = ''
                self.kd      = ''
                self.pn      = ''

        def s(self,string,length):
            string = str(string)
            while len(string)<length: string=' ' + string;
            #if len(string)>length: string=string[0:length];
            return string
        def ll(self,string,length):
            string = str(string)
            while len(string)<length: string=string+' ';
            #if len(string)>length: string=string[0:length];
            return string
        def cprint(self):
            return ' '+self.ll(self.i,4) + self.ll(self.j,4) + self.ll(self.k,4) + self.ll(self.l,4)+s('',2)+s(self.funct,1) + s('%3.1f'%self.phase,10) + s('%2.5f'%self.kd,13) + s((int)(self.pn),6) + '  ; TTT' 

    class impropers:
        def __init__(self,i,j,k,l,funct,phase='',kd='',pn=''):
            self.i       = (i.strip(' '))
            self.j       = (j.strip(' '))
            self.k       = (k.strip(' '))
            self.l       = (l.strip(' '))
            self.funct   = (int)(funct.strip(' '))
            try:
                self.phase   = (float)(phase.strip(' '))
                self.kd      = (float)(kd.strip(' '))
                self.pn      = (float)(pn.strip(' ').strip('\t;').strip('\n'))
            except Exception:
                self.phase   = ''
                self.kd      = ''
                self.pn      = ''

        def s(self,string,length):
            string = str(string)
            while len(string)<length: string=' ' + string;
            #if len(string)>length: string=string[0:length];
            return string
        def ll(self,string,length):
            string = str(string)
            while len(string)<length: string=string+' ';
            #if len(string)>length: string=string[0:length];
            return string

        def cprint(self):
            return self.ll(self.i,4) + self.ll(self.j,4) + self.ll(self.k,4) + self.ll(self.l,4)+s('',5)+s(self.funct,1)+ s('%3.2f'%self.phase,12) + s('%2.5f'%self.kd,12) + s((int)(self.pn),6)+'    ; TTT' 
         

