#!/usr/bin/python
#
#################### W.X. #############################
############### WORMHOLE EXPLORER #####################
#
# 
#
# Wormhole Explorer allows you to study monodromy factorisations
# associated with wormhole pairs.
#
# wx.py introduces two new Python classes:
#
#   HJCF, the class of Hirzebruch-Jung continued fractions
#
#   BOLF, the class of Bhupal-Ozbagci Lefschetz fibrations
#
# Examples of use:
#
# print(BOLF.find(175,64,2))
#
##### prints the (monodromy factorisations for the) Bhupal-Ozbagci
##### Lefschetz fibrations of all Stein fillings of L(175,64) with b_2=2.
#
# print(BOLF.simplify(BOLF.find(175,64))
#
##### prints a simplified version of the monodromy factorisation
##### coming from the (175,64) wormhole
#
#######################################################

from fractions import Fraction
from math import gcd
from pprint import pprint

########## Hirzebruch-Jung continued fractions

class HJCF:
    '''HJCF: The class of Hirzebruch-Jung continued fractions

Usage:

    c=HJCF([c_1,...,c_n])
    
    Returns a HJCF object with Hirzebruch-Jung (HJ) string
    c_1,...,c_n, that is

    c_1-1/(c_2-1/(...-1/c_n)...)

Methods:
    
    c.evaluate()           --- calculate the value of the continued fraction c
    c.find_smaller_zcfs()  --- find zero continued fractions [a_1,...,a_n] with a_i\leq c_i for all i
    c.is_zcf()             --- returns True if c evaluates to zero
    c.is_minimal()         --- returns True if c has no entries equal to 1
    c.blowdown(i)          --- blowdown a continued fraction at position i
    c.find_blowdown_sequence()    --- find a blowdown sequence to a minimal continued fraction (or to [1,1])

Static methods:

    HJCF.calculate_hjcf(a,b)  --- calculate the continued fraction of a/b
    HJCF.hjcf(a,b)            --- equivalent to HJCF(HJCF.calculate_hjcf(a,b))
    '''
    def __init__(self,hjstring):
        '''Usage:

        c=HJCF([c_1,...,c_n])

        Returns a HJCF object with Hirzebruch-Jung (HJ) string
        [c_1,...,c_n], , that is

        c_1-1/(c_2-1/(...-1/c_n)...)

        '''
        self.hjstring=hjstring

    def __repr__(self):
        'Prints a HJCF object'
        return str(self.hjstring)
        
    def evaluate(self):
        '''Usage:

        c.evaluate()

        Computes c_1-1/(c_2-1/...) where c_i are the coefficients of the
        HJ continued fraction.

        '''
        c=self.hjstring
        ans=c[len(c)-1]
        for i in reversed(c[0:len(c)-1]):
            if ans==0:
                return None
            else:
                ans=i-Fraction(1,ans)
        return ans

    def find_smaller_zcfs(self,l):
        '''Usage:

        c.find_smaller_zcfs(l)

        Returns a list of pairs (d,I) where

        - d is a zero-continued fraction, generated from [1,1] by
          blowing up.
        
        - I is a list of length l comprising indices (possibly
          repeated) of coefficients in d,

        such that if, for each i in I, you augment d_i by 1 then you
        get c (so if i is repeated k times then d_i is augmented by k).

        If c=P/Q then these triples correspond under Lisca's
        classification to Stein fillings of the lens space L(P,P-Q).

        '''
        c=self.hjstring # c=[c_1,...,c_n]
        
        # Define an /i-decremented continued fraction/ < c to be a
        # pair [z,(k_1,...,k_i)] where z is obtained from c by
        # decreasing z[k_m] by 1 for each m between 1 and i. We
        # further require that the sequence k_1,...,k_i is decreasing
        # (non-strictly).

        # We generate a tuple of all l-decremented continued fractions
        # < c by an iterative process:

        zlist=([c,()],)
        for i in range(0,l):
            newlist=tuple()
            for a in zlist:
                z=a[0]
                idces=a[1]
                if idces==():
                    bound=len(c)-1
                else:
                    bound=idces[-1]
                for j in range(0,bound+1):
                    if z[j]>=2:
                        decrement=tuple(z[k] if k!=j else z[k]-1 for k in range(0,len(z)))
                        newelt=(decrement,idces+(j,))
                        newlist=newlist+(newelt,)
            zlist=newlist

        # Finally, we test which of the l-decremented continued
        # fractions < c is a zero continued fraction and add those to
        # our list:
        
        smaller_zcfs=[]
        for a in zlist:
            z=HJCF(list(a[0]))
            minimal_model=z.find_blowdown_sequence()[0].hjstring
            if minimal_model==[1,1]:
                smaller_zcfs.append((HJCF(list(a[0])),a[1]))
        return smaller_zcfs

    def blowdown(self,i):
        '''Usage:
        
        c.blowdown(i)

        Given a HJCF object with HJ string c=[c_1,...,c_n] such that
        c_i=1, this function "blows down" at position i, returning the
        new continued fraction [c_1,..,c_{i-1}-1,c_{i+1}-1,...,c_n]

        '''
        c=self.hjstring

        if c==[1]:
            return HJCF([])
        if i==0:
            newhj=[c[1]-1]+c[2:]
        elif i==len(c)-1:
            newhj=c[0:i-1]+[c[i-1]-1]
        elif i==len(c)-2:
            newhj=c[0:i-1]+[c[i-1]-1,c[i+1]-1]
        else:
            newhj=c[0:i-1]+[c[i-1]-1,c[i+1]-1]+c[i+2:]

        return HJCF(newhj)

    def is_minimal(self):
        '''Usage:

        c.is_minimal()

        Returns False if c has a coefficient equal to 1 and False
        otherwise.

        '''
        if not [val for val in self.hjstring if val==1]:
            return True
        else:
            return False
    
    def find_blowdown_sequence(self):
        '''Usage:

        c.find_blowdown_sequence()

        Returns a minimal model for c together with a sequence
        [a_1,...,a_k] such that blowing down c at position a_k,
        a_{k-1},...,a_1 results in either:

        - a minimal continued fraction (if c.evaluate()!=0)

        - or [1,1] (if c.evaluate()==0)

        The algorithm proceeds by blowing down the left-most 1 in c
        and storing its position in a list, then blowing down the
        left-most 1 in the result and appending this on the right etc
        until you hit a minimal guy. There may be other blow-down
        sequences!

        '''
        if self.is_minimal() or self.hjstring==[1,1]:
            return (self,[])
        else:
            # The next line locates the first instance j of the value 1 in
            # the HJ string:
            j,k=next(((i,val) for i,val in enumerate(self.hjstring) if val==1))
            newhjcf=self.blowdown(j)
            nextterm=newhjcf.find_blowdown_sequence()
            return (nextterm[0],nextterm[1]+[j])


    @staticmethod
    def calculate_hjcf(a,b):
        '''Usage:
        
        HJCF.calculate_hjcf(a,b)

        Computes the HJ continued fraction of a/b as a list of
        coefficients.

        '''
        c=[]
        if a%b==0:
            return [a//b]
        else:
            hj=a//b+1
            nxt=hj-Fraction(a,b)
            return [hj]+HJCF.calculate_hjcf(nxt.denominator,nxt.numerator)

    @staticmethod
    def hjcf(a,b):
        '''Usage:

        HJCF.hjcf(a,b)

        Returns a HJCF object representing the continued fraction of
        a/b. Equivalent to HJCF.HJCF(calculate_hjcf(a,b)).

        '''
        return HJCF(HJCF.calculate_hjcf(a,b))

################ Bhupal-Ozbagci Lefschetz Fibrations (BOLFs)
    
class BOLF:
    '''BOLF.

    The class of Bhupal-Ozbagci Lefschetz fibrations (BOLFS).

Background:

    Let D be the closed disc of radius 2 and let D_i be a small open
    disc centred at e^{2\pi i/n}. Let F_n = D\(D_0 u ... u D_{n-1}).

    Given a subset A of {0,1,...,n-1}, there is a unique isotopy class
    of convex simple closed curves C_A in F_n such that C_A encircles
    D_i if and only if i \in A. These are called /convex loops/.

    A Lefschetz fibration over the disc is called /planar/ if its page
    is F_n for some n.

    If we additionally pick a set of vanishing paths and write
    V_1,...,V_k for the corresponding vanishing cycles then we say
    that the Lefschetz fibration has /convex monodromy/ if each V_i is
    a convex loop.
    
    Given:

    - a lens space L(a,b) such that a/(a-b)=[c_1,..,c_n]

    - a zero continued fraction z < c

    there is a Stein filling of L(a,b) (with its Milnor fillable
    contact structure) which comes from the Milnor fibre of a
    Q-Gorenstein smoothing of a P-resolution of the cyclic quotient
    singularity 1/a(1,b).

    Bhupal and Ozbacgi (2016) gave an algorithm for writing down a planar
    Lefschetz fibration with convex monodromies such that the total
    space is this Stein filling.

References:
    
    M. Bhupal and B. Ozbagci, Symplectic fillings of lens spaces as
    Lefschetz fibrations, Journal of the European Mathematical
    Society, Volume 18, Number 7 (2016), 1515--1535
    (or arXiv:1307.6935)

Usage:
        
    F=BOLF(punctures,word,zcf,augmentations)

    punctures --- a list [0,...,n-1] so that the planar fibre of
        our Lefschetz fibration is D - (D_0 u ... u D_{n-1})

    word --- an ordered list of sublists of punctures; each sublist of
        punctures is thought of as specifying a convex loop in the
        fibre and these loops are ordered according to the
        anticlockwise ordering of vanishing cycles for a suitable
        choice of vanishing paths in the base of the Lefschetz
        fibration.

    zcf, augmentations --- the data needed by the Bhupal-Ozbagci
        algorithm to produce a Lefschetz fibration.

Methods:

    F.blowup(i)  --- produces new Lefschetz fibration (corresponds to
        blowing up the zero continued fraction).

    F.augment(i) --- adds a twist to monodromy factorisation.

    F.blowup_sequence()  --- generates a Lefschetz fibration associated to a
        zero continued fraction by a specified sequence of blow-ups.

Static methods:

    BOLF.find(a,b,l) --- generates the BOLFs for all zero continued
        fractions [a_1,...,a_n] with a_i\leq c_i and \sum(c_i-a_i)=l,
        with [c_1,...,c_n]=a/(a-b).

    BOLF.initial() --- generates the Lefschetz fibration for
        the zero continued fraction [1,1]

    '''
    def __init__(self,punctures,word,zcf,augmentations):
        '''Usage:
        
        F=BOLF(punctures,word,zcf,aumentations)

        punctures --- a list [0,...,n-1] so that the planar fibre of
            our Lefschetz fibration is D - (D_0 u ... u D_{n-1})
        
        word --- an ordered list of sublists of punctures; each
            sublist of punctures is thought of as specifying a convex
            loop in the fibre and these loops are ordered according to
            the anticlockwise ordering of vanishing cycles for a
            suitable choice of vanishing paths in the base of the
            Lefschetz fibration.

        zcf, augmentations --- the data needed by the Bhupal-Ozbagci
            algorithm to produce a Lefschetz fibration.

        '''
        self.punctures=punctures
        self.word=word
        self.zcf=zcf
        self.augmentations=augmentations

    def __repr__(self):
        '''Prints information about Lefschetz fibrations.'''
        zcf=[str(z) for z in self.zcf]
        for i in self.augmentations:
            zcf[i]=zcf[i]+"*"
        string="\nZero-continued fraction: |"
        for z in zcf[:-1]:
            string=string+z+","
        string=string+zcf[-1]+"|\n"
        string=string+"Monodromy factorisation:\n"+str(self.word)+"\n\n"
        return string

    def is_augmented(self):
        '''Usage:

        F.is_augmented()

        Returns True if F has been augmented and False otherwise.
        '''
        if self.augmentations:
            return True
        else:
            return False
        
    def blowup(self,i):
        '''Usage:
        
        F.blowup(i)

        If F is a Lefschetz fibration with puncture set [0,...,n-1]
        and monodromy factorisation given by a list of convex Dehn
        twists then F.blowup(i) is a new Lefschetz fibration:
        
        - the ith puncture "buds" into two; all listed convex twists
          which formerly contained i now contain both i and i+1 (and
          the indices from i+1 onwards are shifted up by 1).

        - an extra twist in a new convex loop is added; this loop
          encloses [0,...,i-1,i+1]

        '''
        if self.is_augmented():
            raise TypeError("Cannot blowup an augmented BOLF")

        # There is an additional puncture:
        
        newpunctures=[i for i in range(0,len(self.punctures)+1)]

        # Extra twist in the convex loop enclosing [0,...,i-1,i+1]
        
        extratwist=[k for k in range(0,i)]+[i+1]

        # New monodromy factorisation:
        
        # The function bud tells each convex loop from the previous
        # monodromy factorisation whether it should include i and i+1
        # or not. The new factorisation is obtained from the previous
        # one by budding and adding the extra twist.
        
        bud = lambda x, j: [k for k in x if k<=j]+[k+1 for k in x if
        k>=j]        
        newword=[bud(x,i) for x in self.word]+[extratwist]

        # New zero continued fraction and augmentation
        
        newzcf=self.zcf[0:i-1]+[self.zcf[i-1]+1,1,self.zcf[i]+1]+self.zcf[i+1:len(self.zcf)]
        newaug=()
        
        return BOLF(newpunctures,newword,newzcf,newaug)

    def augment(self,i):
        '''Usage:

        F.augment(i)

        Returns a new Lefschetz fibration with a extra twist in the
        convex loop enclosing [0,...,i].
        
        '''
        extratwist=[j for j in range(0,i+1)]
        newaug=self.augmentations+(i,)
        return BOLF(self.punctures,self.word+[extratwist],self.zcf,newaug)
    
    def blowup_sequence(F,seq):
        '''Usage:
        
        F.blowup_sequence(seq)

        Returns a Lefschetz fibration obtained from F by a specified
        sequence of blow ups.

        '''
        G=F
        for i in seq:
            G=G.blowup(i)
        return G

    @staticmethod
    def find(a,b,l):
        '''Usage:

        BOLF.find(a,b,l)

        Returns the Bhupal-Ozbagci Lefschetz fibrations for any Stein
        filling of L(a,b) with second Betti number l.

        '''
        x=HJCF.hjcf(a,a-b)
        z=x.find_smaller_zcfs(l)
        fibrations=[]
        for zcf in z:
            c,I=zcf
            minimal_model,seq=c.find_blowdown_sequence()
            F=BOLF.initial().blowup_sequence(seq)
            for i in I:
                F=F.augment(i)
            fibrations.append(F)
        return fibrations
    
    @staticmethod
    def mergable(i,j,L):
        '''Here, L is a list of sublists of {0,...,n-1} and we say that the
        triple (i,j,L) is mergable if, for every sublist T in L,
        whenever an element from the subrange i...j appears in T, the
        whole subrange appears in T.

        '''
        R=[x for x in range(i,j+1)]
        for T in L:
            test=[x for x in R if x in T]
            if test!=R and test:
                return False
        return True

    @staticmethod
    def merge(T,i,j):
        return [x for x in T if x<=i]+[x-j+i for x in T if x>j]
    
    @staticmethod
    def iterated_merge(M,N,n):
        '''Usage:
        
        BOLF.iterated_merge(M,N,n)
        
        Given a pair of BOLF monodromy words for the surface F_n, this
        algorithm:

        - strips off any coincident parts of head of tail of the monodromy words,

        - checks if any punctures can be merged without crossing vanishing cycles.

        '''
        while M[0]==N[0]:
            M.pop(0)
            N.pop(0)
        while M[-1]==N[-1]:
            M.pop()
            N.pop()
        mergability=[False,0,0]
        for i in range(0,n):
            for j in reversed(range(i+1,n)):
                if BOLF.mergable(i,j,M+N):
                    mergability=[True,i,j]
                    break
            if mergability[0]:
                break
        if mergability[0]:
            i,j=mergability[1:]
            k=j-i
            K=[BOLF.merge(T,i,j) for T in M]
            L=[BOLF.merge(T,i,j) for T in N]
            return BOLF.iterated_merge(K,L,n-k)
        return M,N 

    @staticmethod
    def simplify(bolfs):
        '''Usage:

        BOLF.simplify(bolfs)

        Here, bolfs is a tuple comprising two BOLFs which are Stein
        fillings of the same lens space, this returns a simplification
        of the monodromy substitution to go from F to G.

        '''
        F,G=bolfs
        if F==G:
            raise ValueError("These are the same BOLF.")
        return BOLF.iterated_merge(F.word,G.word,len(F.punctures))
    
    @staticmethod
    def initial():
        '''Usage:

        F=BOLF.initial()

        Creates a Lefschetz fibration whose page is a twice punctured
        disc and whose monodromy is one twist around the second
        puncture. This corresponds to the zero continued fraction
        [1,1]; all other Lefschetz fibrations corresponding to zero
        continued fractions can be obtained by a sequence of blow-ups
        from this one.

        '''
        return BOLF([0,1],[[1]],[1,1],())

def find_wormholes(N):
    '''Usage:

    find_wormholes(N)

    Returns a dictionary whose keys are all lens spaces L(p,q) with
    p<N admitting two Stein fillings with b_2=1, and whose entry for
    the key L(p,q) is the pair of augmented zero continued fractions
    giving the fillings.

    '''
    wormholes={}
    for p in range(2,N):
        qs=[q for q in range(1,p//2+1) if gcd(p,q)==1]
        for q in qs:
            L=HJCF.hjcf(p,p-q).find_smaller_zcfs(2)
            if len(L)==2:
                wormholes[(p,q)]=L
    wormholes_strings={}
    for keys in wormholes:
        J,K=wormholes[keys]
        Jzcf=[str(z) for z in J[0].hjstring]
        for i in J[1]:
            Jzcf[i]=Jzcf[i]+"*"
        Kzcf=[str(z) for z in K[0].hjstring]
        for i in K[1]:
            Kzcf[i]=Kzcf[i]+"*"
        wormholes_strings[keys]=(Jzcf,Kzcf)
    return wormholes_strings

def print_wormholes(N):
    '''Usage:

    print_wormholes(N)

    Prints the output of find_wormholes(N) in a more readable way.

    '''
    wm=find_wormholes(N)
    for keys in wm:
        print(keys,":")
        for filling in wm[keys]:
            print('   [%s]' % ', '.join(map(str, filling)))

