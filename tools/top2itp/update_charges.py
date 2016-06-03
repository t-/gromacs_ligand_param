from gmx_top import *
import os
import sys

newtop=sys.argv[1]
oldtop=sys.argv[2]

ntop =topology()
otop =topology()                    
ntop.read_top(newtop)
otop.read_top(oldtop)
print 'new',len(ntop.atoms_list)
print 'old',len(otop.atoms_list)

for i in xrange(len(ntop.atoms_list)):
    otop.atoms_list[i].charge=ntop.atoms_list[i].charge

otop.write_top(oldtop)
