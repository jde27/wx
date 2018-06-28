#!/usr/bin/python

from wx import *
from epr import *
from operator import itemgetter
from itertools import groupby

lens_spaces=epr_by_delta(30,30,10)
monod=[]

f = lambda a,b: BOLF.simplify(BOLF.find(a,b,2))

for lens in lens_spaces[2]:
    a,b=lens
    monod.append({'lens':lens,'factor':f(a,b)})
monod.sort(key=itemgetter('factor'))
for factor,items in groupby(monod,key=itemgetter('factor')):
    print(factor)
    for i in items:
        print('    ',i['lens'])
