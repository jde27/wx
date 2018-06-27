#!/usr/bin/python
#
############### WORMHOLE EXPLORER #####################
#
# Wormhole Explorer allows you to study monodromy factorisations
# associated with wormhole pairs.
#
# wexp.py introduces two new Python classes:
#
# CF, the class of continued fractions
#
# Fibration, the class of (planar) Lefschetz fibrations whose
# monodromies are convex Dehn twists
#
#######################################################

from fractions import Fraction

class CF:
    '''CF: The class of continued fractions

Usage:

    c=CF([c_1,...,c_n])
    
    Returns a CF object with Hirzebruch-Jung (HJ) string
    c_1,...,c_n.

Methods:
    
    c.evaluate()   --- calculate the value of the continued fraction c
    c.find_zcfs()  --- find zero continued fractions [a_1,...,a_n] with a_i\leq c_i for all i
    c.blowdown(i)  --- blowdown a zero continued fraction
    c.findseq()    --- find a blowdown sequence for a zero-continued fraction

Static methods:

    CF.calculate_cf(a,b)  --- calculate the continued fraction of a/b
    CF.cf(a,b)            --- equivalent to CF(CF.calculate_cf(a,b))
    '''
    def __init__(self,hjstring):
        '''Usage:

        c=CF([c_1,...,c_n])

        Returns a CF object with Hirzebruch-Jung (HJ) string
        c_1,...,c_n.

        '''
        self.hjstring=hjstring

    def evaluate(self):
        '''Usage:

        c.evaluate()

        Computes c_1-1/(c_2-1/...) where c_i are the coefficients of the
        continued fraction.

        '''
        c=self.hjstring
        ans=c[len(c)-1]
        for i in reversed(c[0:len(c)-1]):
            ans=i-Fraction(1,ans)
        return ans

    def find_zcfs(self):
        '''Usage:

        c.find_zcfs()

        Returns a list of 3-tuples (d,i,j) where

        - d is a zero-continued fraction,
        
        - i and j are indices of coefficients in d

        such that if you augment d_i and d_j each by 1 then you get c

        If c=P/Q then these triples correspond under Lisca's
        classification to Stein fillings of the lens space L(P,P-Q).

        '''
        c=self.hjstring # c=[c_1,...,c_n]
        idx=[]
        for i in range(0,len(c)):
            for j in range(i+1,len(c)):
                # d = [c_1,...,c_i-1,...,c_j-1,...,c_n]
                d=CF([c[k] if k not in [i,j] else c[k]-1 for k in range(0,len(c))])                
                if d.evaluate()==0:
                    idx.append((d,i,j)) # if d is a zero continued
                                        # fraction, add this 3-tuple
                                        # to the list
        return idx

    def blowdown(self):
        '''Usage:
        
        c.blowdown()

        Given a continued fraction c=[c_1,...,c_n], this function
        finds the smallest j such that c_j=1 and "blows down" at that
        point, replacing it with a continued fraction
        [c_1,...,c_{j-1}-1,c_{j+1}-1,...c_n]. If the continued
        fraction does not contain a 1, this raises a type error.

        '''
        s=self.hjstring
        # The next line locates the first instance of the value 1 in
        # the HJ string
        j,k=next(((i,val) for i,val in enumerate(s) if val==1),(-1,1))
        if j==-1:
            raise TypeError("This continued fraction is minimal and cannot be blown down.")
        if j==0: # This if..elif sequence performs the blowdown
                 # depending on whereabouts you are blowing down
            newhj=[s[1]-1]+s[2:]
        elif j==len(s)-1:
            newhj=s[0:j-1]+[s[j-1]-1]
        elif j==len(s)-2:
            newhj=s[0:j-1]+[s[j-1]-1,s[j+1]-1]
        else:
            newhj=s[0:j-1]+[s[j-1]-1,s[j+1]-1]+s[j+2:]
        return j,newhj 

    def findseq(self):
        '''Usage:

        c.findseq()

        Returns a sequence [a_1,a_2,...,a_k] such that starting from
        the zero-continued fraction [1,1] and blowing up at position
        a_1, then position a_2, etc results in the zero-continued
        fraction c. If c is not a zero continued fraction then this
        raises a type error.

        '''
        if self.evaluate()==0: # Check it's a zero-continued fraction
            if self.hjstring==[1,1]: # If we've reached [1,1], there's
                                     # nothing more to do
                return []
            else:
                i,newhj=self.blowdown() # otherwise, blow down and
                                        # record the position you blew
                                        # down...
                return CF(newhj).findseq()+[i] # ...and continue
        else:
            raise TypeError("Cannot find blowdown sequence for a nonzero continued fraction.")

    @staticmethod
    def calculate_cf(a,b):
        '''Usage:
        
        CF.calculate_cf(a,b)

        Computes the HJ continued fraction of a/b as a list of
        coefficients.

        '''
        c=[]
        if a%b==0:
            return [a//b]
        else:
            hj=a//b+1
            nxt=hj-Fraction(a,b)
            return [hj]+CF.calculate_cf(nxt.denominator,nxt.numerator)

    @staticmethod
    def cf(a,b):
        '''Returns a CF object representing the continued fraction of a/b.

        '''
        return CF(CF.calculate_cf(a,b))

                    
class Fibration:
    '''Fibration.

    The class of planar Lefschetz fibrations whose vanishing cycles are
    convex loops.

    We consider the fibre to be the closed disc D of radius 2 minus
    small open discs D_0,...,D_{n-1} centred at the nth roots of
    unity. A loop in D - (D_0 u ... u D_{n-1}) is called /convex/ if
    it is the boundary of a convex region M in D; we can specify a
    unique isotopy class of convex loops by specifying which D_i are
    contained in the interior of M, so we can specify the vanishing
    cycles of our Lefschetz fibration just by giving a subset of
    [0,...,n-1].

Usage:
        
    F=Fibration(punctures,word,zcf)

    punctures --- a list [0,...,n-1] so that the planar fibre of
        our Lefschetz fibration is D - (D_0 u ... u D_{n-1})

    word --- an ordered list of sublists of punctures; each sublist of
        punctures is thought of as specifying a convex loop in the
        fibre and these loops are ordered according to the
        anticlockwise ordering of vanishing cycles for a suitable
        choice of vanishing paths in the base of the Lefschetz
        fibration.

    zcf --- the zero continued fraction corresponding to this
        Lefschetz fibration. # currently this is here, but I may
        remove it later, as there are planar Lefschetz fibrations that
        don't come from zero continued fractions.

Methods:

    F.blowup(i)  --- produces new Lefschetz fibration (corresponds to
        blowing up the zero continued fraction).

    F.augment(i) --- adds a twist to monodromy factorisation.

    F.blowseq()  --- generates a Lefschetz fibration associated to a
        zero continued fraction by a specified sequence of blow-ups.

Static methods:

    find(a,b) --- generates the Lefschetz fibrations for all zero
        continued fractions [a_1,...,a_n] with a_i\leq c_i and
        [c_1,...,c_n]=a/(a-b).

    Fibration.initial() --- generates the Lefschetz fibration for
        the zero continued fraction [1,1]

    '''
    def __init__(self,punctures,word,zcf):
        '''Usage:
        
        F=Fibration(punctures,word,zcf)

        punctures --- a list [0,...,n-1] so that the planar fibre of
        our Lefschetz fibration is D - (D_0 u ... u D_{n-1})

        word --- an ordered list of sublists of punctures; each
        sublist of punctures is thought of as specifying a convex loop
        in the fibre and these loops are ordered according to the
        anticlockwise ordering of vanishing cycles for a suitable
        choice of vanishing paths in the base of the Lefschetz
        fibration.

        zcf --- the zero continued fraction corresponding to this
        Lefschetz fibration. # currently this is here, but I may
        remove it later, as there are planar Lefschetz fibrations that
        don't come from zero continued fractions.

        '''
        self.punctures=punctures
        self.word=word
        self.zcf=zcf

    def __repr__(self):
        '''For printing information about Lefschetz fibrations.'''
        string="Zero-continued fraction "+str(self.zcf)+": \nMonodromy factorisation \n"+str(self.word)
        return string
        
    def blowup(self,i):
        '''Usage:
        
        F.blowup(i)

        If F is a Lefschetz fibration with puncture set [0,...,n-1]
        and monodromy factorisation given by a list of convex Dehn
        twists then F.blowup(i) is a new Lefschetz fibration:
        
        - the ith puncture "buds" into two; all listed convex twists
        which formerly contained i now contain both i and i+1 (and the
        indices from i+1 onwards are shifted up by 1).

        - an extra twist in a new convex loop is added; this loop
        encloses [0,...,i-1,i+1]

        '''
        newpunctures=[i for i in range(0,len(self.punctures)+1)] # There is an additional puncture
        bud = lambda x, j: [k for k in x if k<=j]+[k+1 for k in x if
        k>=j] # This function tells each convex loop from the previous
        # monodromy factorisation whether it should include i and i+1
        # or not.
        extratwist=[k for k in range(0,i)]+[i+1] # This is the extra
        # twist in the convex loop enclosing [0,...,i-1,i+1]
        newword=[bud(x,i) for x in self.word]+[extratwist] # New monodromy factorisation
        newzcf=self.zcf[0:i-1]+[self.zcf[i-1]+1,1,self.zcf[i]+1]+self.zcf[i+1:len(self.zcf)]
        # again, I may remove the zcf variable later
        return Fibration(newpunctures,newword,newzcf)

    def augment(self,i):
        '''Usage:

        F.augment(i)

        Returns a new Lefschetz fibration with a extra twist in the
        convex loop enclosing [0,...,i].
        
        '''
        extratwist=[j for j in range(0,i+1)]
        newzcf=["%s*" % self.zcf[j] if j==i else self.zcf[j] for j in
        range(0,len(self.zcf))]
        # again, I may remove the zcf variable later
        return Fibration(self.punctures,self.word+[extratwist],newzcf)
    
    def blowseq(F,seq):
        '''Usage:
        
        F.blowseq(seq)

        Returns a Lefschetz fibration obtained from F by a specified
        sequence of blow ups.

        '''
        G=F
        for i in seq:
            G=G.blowup(i)
        return G

    @staticmethod
    def find(a,b):
        '''Usage:

        Fibration.find(a,b)

        Returns the Lefschetz fibrations for any wormhole pair
        associated to the lens space L(a,b).
        '''
        x=CF.cf(a,a-b)
        z=x.find_zcfs()
        fibrations=[]
        for triplet in z:
            c,i,j=triplet
            seq=c.findseq()
            F=Fibration.initial().blowseq(seq).augment(i).augment(j)
            fibrations.append(F)
        return fibrations

    @staticmethod
    def initial():
        '''Usage:

        F=Fibration.initial()

        Creates a Lefschetz fibration whose page is a twice punctured
        disc and whose monodromy is one twist around the second
        puncture. This corresponds to the zero continued fraction
        [1,1]; all other Lefschetz fibrations corresponding to zero
        continued fractions can be obtained by a sequence of blow-ups
        from this one.

        '''
        return Fibration([0,1],[[1]],[1,1])

