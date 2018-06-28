#!/usr/bin/python

from collections import defaultdict
from pprint import pprint
import itertools
import math

def possible_qs(p,idx):
    '''Returns a list of possible q values given p.

    If p=1 then q can be 0 or 1 (if idx=1 or 2 respectively).

    Otherwise, q can be anything between 1 and p-1 coprime to p.

    '''
    if p==1:
        return [idx-1]
    else:
        return [q for q in range(1,p//2+1) if math.gcd(p,q)==1]

def sdo(T):
    '''Computes Delta, Omega and sigma for the given Wahl singularities
    and shear invariant.

    '''
    (p_1,q_1,p_2,q_2,c)=T
    sigma = (c-1)*p_1*p_2+p_2*q_1-p_1*q_2
    Delta = p_1*p_1+p_2*p_2+sigma*p_1*p_2
    OmegaFull = p_1*q_1+p_2*q_2-1+sigma*p_2*q_1-(c-1)*p_2*p_2
    Omega = OmegaFull % Delta
    return [sigma,Delta,Omega]
    
def enum_epr(P_1,P_2,C):
   '''Enumerates all extremal P-resolutions:

      Enumerates all:
   
        p_1\leq P_1, 0\leq q_1<p_1,
        p_2\leq P_2, 0<q_2\leq p_2,
        c\leq C

      such that the polygon \Pi(p_1,q_1,p_2,q_1,c,*) is a K-positive
      truncation of a wedge \pi(\Delta,\Omega) with
      0<\Omega<\Delta. Recall that K-positivity is defined by:
  
      \sigma(\Pi)=(c-1)p_1p_2+p_2q_1-p_1q_2 > 0

      (we always call |\sigma(\Pi)|=\delta). Note that

      \Delta(\Pi) = p_1^2+p_2^2+\sigma(\Pi)p_1p_2
   
      and

      \Omega(\Pi) = p_1q_1+p_2q_2-1+\sigma(\Pi)p_2q_1+(c-1)p_2^2
                     ( mod \Delta(\Pi) )
   
      The data is stored as a defaultdict (i.e. a dictionary indexed
      by (\Delta,\Omega) whose values form a set of tuples
      (p_1,q_1,p_2,q_2,c)).

   '''
   epr = defaultdict(set)
   for (p_1,p_2) in itertools.product(range(1,P_1+1),range(1,P_2+1)):
           Q_1=possible_qs(p_1,1)
           Q_2=possible_qs(p_2,2)
           for (q_1,q_2) in itertools.product(Q_1,Q_2):
               for c in range (1,C+1):
                   [sigma,Delta,Omega]=sdo((p_1,q_1,p_2,q_2,c))
                   if sigma>0 and Delta>0:
                       epr[(Delta,Omega)].add((p_1,q_1,p_2,q_2,c))
   return epr

def palindrome(T):
    (p_1,q_1,p_2,q_2,c)=T
    if p_1==1:
        (r_2,s_2)=(1,1)
    else:
        (r_2,s_2)=(p_1,p_1-q_1)
    if p_2==1:
        (r_1,s_1)=(1,0)
    else:
        (r_1,s_1)=(p_2,p_2-q_2)
    return (r_1,s_1,r_2,s_2,c)

def epr_by_delta(n1,n2,cbound):
    epr=enum_epr(n1,n2,cbound)
    grp_by_delta=defaultdict(list)
    for keys in epr:
        values=list(epr[keys])
        if len(values)==2:
            delta=sdo(values[0])[0]
            grp_by_delta[delta].append(keys)
    for delta in grp_by_delta:
        grp_by_delta[delta]=sorted(grp_by_delta[delta])
    return grp_by_delta
