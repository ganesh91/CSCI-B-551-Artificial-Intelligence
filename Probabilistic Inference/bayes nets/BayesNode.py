"""
BayesNode represent the basic modular node of the BayesNet.
BayesNode defines the properties of the Random Variable like
indegree, Parent Random Variable, conditional probability table.
"""
from collections import defaultdict
import pprint
import copy

pp = pprint.PrettyPrinter(indent=4)

class BayesNode:
    def __init__(self,id,parents,indegree):
        """
        Initialize a node in Bayes Net.
        """
        self.id=id
        self.parents=None if id==parents else parents.split(",")
        self.indegree=0 if self.parents is None else len(self.parents)
        self.cpt=None
        self.isBoolean=True

    def buildCPT(self):
        """
        Build the conditional probability table for the node.
        Also calculates the negation of the node.
        Saves CPT as dict of dicts.
        Result: for a p(A|B):

        (B) {T: {A:0.00,~A:0.00}
        (B) F: {A:0.00,~A:0.00}}
        """
        if self.indegree==0:
            self.cpt={self.id:None,"~"+self.id:None}
        else:
            queue=copy.copy(self.parents)
            queue.append(self.id)
            branch={'T':None,'F':None}

            def buildCPTRecurse(parentstack,existingstack,ret=None):
                if ret!=None:
                    return existingstack
                if len(parentstack)!=0:
                    parentstack.pop()
                    if len(existingstack)==0:
                        exstk=[]
                        for i in range(2*len(parentstack)):
                            exstk.append({self.id:None,"~"+self.id:None})
                        buildCPTRecurse(parentstack,exstk)
                    else:
                        exstk=[]
                        i=0
                        for i in range(0,len(existingstack),2):
                            branch={'T':existingstack[i],'F':existingstack[i+1]}
                            exstk.append(branch)
                        buildCPTRecurse(parentstack,exstk)
                else:
                    self.cpt=existingstack[0]
                    return existingstack
            buildCPTRecurse(queue,[])

    def setValue(self,compressed,value):
        """
        Parameters: Compressed(concatenated) CPT entry : Eg "T,F"
        value: float/integer
        returns None.
        """
        values=compressed.split(",")
        if self.id==compressed:
            self.cpt={self.id:float(value),"~"+self.id:1-float(value)}
        else:
            initial=True
            while values:
                a=values.pop(0)
                if initial:
                    holder=self.cpt[a]
                    initial=False
                else:
                    holder=holder[a]
            holder[self.id]=float(value)
            holder["~"+self.id]=1-float(value)

    def getValue(self,X,givenY):
        """
        Given the query (Eg, for Node A, valid queries are [A,~A]) append
        evidence (givenY) in compressed form ("T,F,T"), return cpt entry
        """
        givenList=givenY.split(",")
        cpt=self.cpt
        if self.indegree==0:
            if givenY == 'T':
                return(cpt[self.id])
            else:
                return(cpt["~"+self.id])
        else:
            for i in givenList:
                cpt=cpt[i]
            return(cpt[X])

    def print(self):
        """
        print the contents of the Bayes Node.
        """
        #print(self.id,self.indegree,self.outdegree,self.parents,self.children,self.cptenties,self.cpt)
        print("ID :",self.id, "Indegree :",self.indegree, "Parents :",self.parents)
        pp.pprint(self.cpt)
