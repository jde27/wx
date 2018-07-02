#!/usr/bin/python

from wx import *
from epr import *
from operator import itemgetter
from itertools import groupby

delta=2
P1=30
P2=30
C=10

lens_spaces=epr_by_delta(P1,P2,C)
monod=[]

f = lambda a,b: BOLF.simplify(BOLF.find(a,b,2))

for lens in lens_spaces[delta]:
    a,b=lens
    monod.append({'lens':lens,'factor':f(a,b)})
monod.sort(key=itemgetter('factor'))
for factor,items in groupby(monod,key=itemgetter('factor')):
    print(factor)
    for i in items:
        print('    L(',i['lens'],')')
