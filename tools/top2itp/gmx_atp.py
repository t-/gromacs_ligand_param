class gmx_atp():
    def __init__(self,top):
        self.top = top
    def write_atp(self,filename):
        def r(string,length):
            string = str(string)
            while len(string)<length: string=string+' ';
            return string
        def l(string,length):
            string = str(string)
            while len(string)<length: string=' '+string;
            return string
        print 'writing',filename
        f = open(filename,'w')
        atp = {}
        for i in self.top.atoms_list:
            atp[r(i.type,10)+l('%2.5f' %i.mass,16)+str('\t; TTT')]=''
        for i in atp.keys():
            if i.upper()!=i and i.find('amber')<0:
                #if i not in  ['CT','H1','C','O','N','N3','HC'] and i.find('amber')<0:   #only print those entries which are not already known by ffamber99sb --> i.e. do not contain amber99_xx
                print >>f, i
        f.close()
