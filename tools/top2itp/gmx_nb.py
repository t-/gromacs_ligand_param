class ffnb:
    def __init__(self):  
        self.atomtypes      = []
        self.pairtypes          = []
        self.atomtypes_meta     = []
        self.pairtypes_meta         = []
    def write(self,filename):
        print 'writing',filename
        f = open(filename,'w')
        print >>f, '[ atomtypes ]'
        for i in self.atomtypes:
            print >>f, i.cprint()
    def read_top(self,filename):
        block=[]
        blocklist=[]
        #Read Topology and separate it into blocks [ atoms ], [ bonds ], etc.
        for i in open(filename,'r'):
            if len(i.strip('\n'))==0:    # blank line indicates the end of a block [ atoms ], [ bonds ], etc.
                if len(block)>0: blocklist.append(block); 
                block = []
            elif len(i.strip('\n'))>0:    # read block
                block.append(i.strip('\n'))
        blocklist.append(block);
        
        #Read [ atomtypes ]
        read = 'true'
        for block in blocklist:
            if block[0].find('atomtypes') >= 0:
                for i in range(len(block)): #lines
                    if block[i][0]==';' or block[i][0]=='[' or block[i][0].find('#')>=0 :
                        self.atomtypes_meta.append(block[i])           #save comment lines, [...], etc.
                        if block[i].find('#ifdef')>=0:
                            read = 'false'
                        if block[i].find('#endif')>=0:
                            read = 'true'
                    elif read == 'true':
                        top=block[i].split('\t')
                        for i in range(top.count('')): top.remove(''); #remove blanks from array
                        #name,at_num,mass,charge,ptype,sigma,epsilon
                        self.atomtypes.append(self.atomtype(top[0],top[1],top[2],top[3],top[4],top[5],top[6]))
                                            
        
        #Read [ pairs ]
        read = 'true'
        for block in blocklist:
            if block[0].find('pairtypes') >= 0:
                for i in range(len(block)): #lines
                    if block[i][0]==';' or block[i][0]=='[' or block[i][0].find('#')>=0 :
                        self.pairtypes_meta.append(block[i])           #save comment lines, [...], etc.
                        if block[i].find('#ifdef')>=0:
                            read = 'false'
                        if block[i].find('#endif')>=0:
                            read = 'true'
                    elif read == 'true':
                        top=block[i].split('\t')
                        for i in range(top.count('')): top.remove(''); #remove blanks from array
                        #i    j    func    sigma1-4    epsilon1-4
                        self.pairtypes.append(self.pairtype(top[0],top[1],top[2],top[3],top[4]))

    class atomtype:
        def __init__(self,name,atnum,mass,charge,ptype,sigma,epsilon):
            self.name    = name.strip(' ')
            self.atnum  = atnum.strip(' ')
            self.mass    = (float)(mass.strip(' '))
            self.charge  = (float)(charge.strip(' '))
            self.ptype   = ptype.strip(' ')
            self.sigma   = (float)(sigma.strip(' '))
            self.epsilon = (float)(epsilon.split(' ')[0].strip(' ').strip(';')) # kill comments and remove format errors
        def s(self,string,length):
            string = str(string)
            while len(string)<length: string=' ' + string;
            return string
        def l(self,string,length):
            string = str(string)
            while len(string)<length: string= string + ' ';
            return string
        def cprint(self):
            return self.l(str(self.name),10)+self.s(str(self.atnum),4)+self.s(str(self.mass),11)+'  '+self.s('%1.4f'%self.charge,8)+self.s(str(self.ptype),3)+self.s('%e'%self.sigma,14)+self.s('%e'%self.epsilon,13)
    
    class pairtype:
        #not needed so far... if you do, plz contact me
        def __init__(self,ai,aj,funct,sigma14,epsilon14):
            self.ai        = (ai.strip(' '))
            self.aj        = (aj.strip(' '))
            self.funct     = (int)(funct.strip(' '))
            self.sigma14   = (float)(sigma14.strip(' '))
            self.epsilon14 = (float)(epsilon14.strip(' '))
        def cprint(self):
            return 'broken please fix'
            #return str(self.ai)+'\t'+str(self.aj)+'\t'+str(self.funct)+'\t'+str(self.sigma14)+'\t'+str(self.epsilon14)
                        