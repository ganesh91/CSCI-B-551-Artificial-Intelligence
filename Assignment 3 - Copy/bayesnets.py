from collections import defaultdict
from BayesNode import BayesNode
import copy
import pprint
import random
from functools import reduce
import operator

pp = pprint.PrettyPrinter(indent=4)

class BayesNets:
    def __init__(self):
        """
        Class Member Declarations
        """
        self.adjacencyMetadata=defaultdict(str)
        self.adjacency=defaultdict(str)
        self.topology={'B':1,'E':2,'A':3,'J':4,'M':5}

    def returnIndexes(self,X,array):
        if X=='B':
            return array[0]
        if X=='E':
            return array[1]
        if X=='A':
            return ",".join(array[0:2])
        if X=='J':
            return array[2]
        if X=='M':
            return array[2]

    def print(self):
        pp.pprint(self.adjacencyMetadata)
        for item in self.adjacency.keys():
            self.adjacency[item].print()

    def readAdjacencyGraph(self,filename):
        """
        Create Adjacency Matrix representation of the graph
        """
        try:
            for line in open(filename,'r'):
                incoming,outgoing=line.strip().split(":")
                no_outgoing=outgoing.split(",")
                self.adjacencyMetadata[incoming]=dict(zip(no_outgoing,range(len(no_outgoing))))
                if incoming not in self.adjacency.keys():
                    self.adjacency[incoming]=None
                for item in no_outgoing:
                    if item not in self.adjacency.keys():
                        self.adjacency[item]=None
        except Exception as e:
            raise

    def buildAdjacencyMetdatata(self,filename):
        """
        Given the adjacency list, read through the possible values
        create the bayes net.
        """
        initial=True
        node=None
        values=[]
        for line in open(filename,'r'):
            params=line.strip().split(":")
            if len(params)==4:
                randomvariable=params[0]
                parents=params[2]
                tablelen=params[3]
                if initial:
                    initial=False
                    node=BayesNode(params[0],params[2],params[3])
                    node.buildCPT()
                else:
                    self.adjacency[node.id]=node
                    if params[0]!="":
                        node=BayesNode(params[0],params[2],params[3])
                        node.buildCPT()
            else:
                node.setValue(params[0],params[1])

    def inferTopology(self,variables):
        variableList=[(i,self.topology[i]) for i in variables]
        variableList=sorted(variableList,key=lambda x: x[1])
        return ([a for (a,b) in variableList])

    def inferHiddenVariables(self,variables):
        topology=[(item,self.topology[item]) for item in self.topology]
        topology=[a for (a,b) in topology if b<=self.topology[variables[-1]]]
        return (topology)

    def expandNode(self,query,evidence):
        querynodes=set(query)
        evidencenodesgiven=set([a for (a,b) in evidence])
        nodes=querynodes.union(evidencenodesgiven)
        orderednodes=self.inferTopology(self.inferHiddenVariables(self.inferTopology(nodes)))
        evidencenodes=self.inferTopology(self.inferHiddenVariables(self.inferTopology(evidencenodesgiven)))
        return(querynodes,evidencenodesgiven,orderednodes,evidencenodes)

    def querytovector(self,topology,query,evidence,negate=False):
        vector=[]
        for item in topology:
            if item in query:
                if not negate:
                    vector.append('T')
                else:
                    vector.append('F')
            elif item in [a for (a,b) in evidence]:
                vector.append([b for (a,b) in evidence if a==item][0])
            else:
                vector.append('-')
        return vector

    def vectorcount(self,sample,query):
        same=0
        for item in sample:
            local=0
            for i in range(len(item)):
                    if query[i]==item[i]:
                        local+=1
                    elif query[i] == '-':
                        local+=1
            if local==len(query):
                same+=1
        #print(query,same)
        return same

    def gibbsvectorcount(self,sample,query):
        #print(sample[0][0],query)
        same=[]
        for item in sample:
            local=0
            for i in range(len(item[0])):
                    if query[i]==item[0][i]:
                        local+=1
                    elif query[i] == '-':
                        local+=1
            if local==len(query):
                same.append(item[1])
        same=sum(same)
        #print(query,same)
        return same

    def sampledistribution(self,nsample,orderednodes):
        samples=[]
        #print(orderednodes)
        j=0
        while j<nsample:
            local=[]
            for i in orderednodes:
                indegree=self.adjacency[i].indegree
                if indegree==0:
                    cp=self.adjacency[i].getValue(i,'T')
                    ra=random.random()
                    #print(cp,ra)
                    if ra < cp :
                        local.append('T')
                    else:
                        local.append('F')
                else:
                    cp=self.adjacency[i].getValue(i,self.returnIndexes(i,local))
                    ra=random.random()
                    #print(cp,ra)
                    if ra < cp :
                        local.append('T')
                    else:
                        local.append('F')
            samples.append(local)
            #print(local)
            j+=1
        return(samples)

    def priorsampling(self,nsample,query,evidence):
        closedquery=[(qt,qt) for qt in query]
        closedevidence=evidence
        openevidence=[a for (a,b) in evidence]
        querynodes,evidencenodesgiven,orderednodes,evidencenodes=self.expandNode(query,evidence)
        query=self.querytovector(orderednodes,query,evidence)
        evidence=self.querytovector(evidencenodes,query,evidence)
        samples=self.sampledistribution(nsample,orderednodes)
        for _ in range(len(query)-len(evidence)):
            evidence.append("-")
        a=self.vectorcount(samples,query)
        b=self.vectorcount(samples,evidence)
        if a==0 or b==0:
            return(0)
        else:
            return(a/b)

    def rejectsampledistribution(self,nsample,orderednodes,evidence):
        samples=[]
        #print(orderednodes)
        #print(evidence)
        j=0
        while j<nsample:
            local=[]
            count=0
            for num,i in enumerate(orderednodes):
                indegree=self.adjacency[i].indegree
                if indegree==0:
                    cp=self.adjacency[i].getValue(i,'T')
                    ra=random.random()
                    if ra < cp :
                        if evidence[num]=='T' or evidence[num]=='-':
                            local.append('T')
                        else:
                            break
                    else:
                        if evidence[num]=='F' or evidence[num]=='-':
                            local.append('F')
                        else:
                            break
                else:
                    if len(local) <= 4:
                        s=",".join(local[-indegree:])
                    else:
                        s=",".join(local[0:4])
                    cp=self.adjacency[i].getValue(i,self.returnIndexes(i,local))
                    ra=random.random()
                    if ra < cp :
                        if evidence[num]=='T' or evidence[num]=='-':
                            local.append('T')
                        else:
                            break
                    else:
                        if evidence[num]=='F' or evidence[num]=='-':
                            local.append('F')
                        else:
                            break
                if len(local)==len(orderednodes):
                    #print(local,evidence)
                    samples.append(local)
                    j+=1
        return(samples)

    def rejectionsampling(self,nsample,query,evidence):
        closedquery=[(qt,qt) for qt in query]
        closedevidence=evidence
        openevidence=[a for (a,b) in evidence]
        querynodes,evidencenodesgiven,orderednodes,evidencenodes=self.expandNode(query,evidence)
        query=self.querytovector(orderednodes,query,evidence)
        evidence=self.querytovector(evidencenodes,query,evidence)
        for _ in range(len(query)-len(evidence)):
            evidence.append("-")
        samples=self.rejectsampledistribution(nsample,orderednodes,evidence)
        a=self.vectorcount(samples,query)
        b=self.vectorcount(samples,evidence)
        if a==0 or b==0:
            return(0)
        else:
            return(a/b)

    def mldistribution(self,nsample,orderednodes,evidence):
        samples=[]
        #print(orderednodes)
        #print(evidence)
        j=0
        while j<nsample:
            local=[]
            weight=[]
            count=0
            for num,i in enumerate(orderednodes):
                if evidence[num]=='F':
                    qv="~"+i
                else:
                    qv=i
                indegree=self.adjacency[i].indegree
                if indegree==0:
                    if evidence[num]=='-':
                        cp=self.adjacency[i].getValue(qv,'T')
                        ra=random.random()
                        if ra < cp :
                            local.append('T')
                        else:
                            local.append('F')
                    else:
                        weight.append(self.adjacency[i].getValue(qv,evidence[num]))
                        local.append(evidence[num])
                else:
                    if len(local) <= 4:
                        s=",".join(local[-indegree:])
                    else:
                        s=",".join(local[0:4])
                    cp=self.adjacency[i].getValue(i,s)
                    ra=random.random()
                    if evidence[num]=='-':
                        if ra < cp :
                            local.append('T')
                        else:
                            local.append('F')
                    else:
                        weight.append(self.adjacency[i].getValue(qv,self.returnIndexes(i,local)))
                        local.append(evidence[num])
                if len(local)==len(orderednodes):
                    #print(local,evidence,weight)
                    samples.append((local,reduce(operator.mul,weight,1)))
            j+=1
        return(samples)

    def enumeration(self,query,evidence):
        order=self.topology.items()
        order=[a for (a,b) in sorted(order,key=lambda x: x[1])]
        possible=['T','F']
        jointprobs=[]
        for B in possible:
            for E in possible:
                for A in possible:
                    for J in possible:
                        for M in possible:
                            enumerated_list=[B,E,A,J,M]
                            product=1
                            for q,variable in enumerate(order):
                                if enumerated_list[q]=='F':
                                    qv="~"+variable
                                else:
                                    qv=variable
                                product=product*self.adjacency[variable].getValue(qv,self.returnIndexes(variable,enumerated_list))
                            jointprobs.append((enumerated_list,product))
        closedquery=[(qt,qt) for qt in query]
        closedevidence=evidence
        openevidence=[a for (a,b) in evidence]
        querynodes,evidencenodesgiven,orderednodes,evidencenodes=self.expandNode(query,evidence)
        query2=self.querytovector(orderednodes,query,evidence)
        query1=self.querytovector(orderednodes,query,evidence,True)
        #query1=self.querytovector(orderednodes,query,evidence,True)
        #evidence=self.querytovector(orderednodes,query,evidence)
        for _ in range(len(jointprobs[0])-len(query1)):
            query1.append("-")
            query2.append("-")
        #print(query1,query2)
        a=self.gibbsvectorcount(jointprobs,query2)
        b=self.gibbsvectorcount(jointprobs,query1)
        if a==0 or b==0:
            return(0)
        else:
            return(a/(a+b))

    def maxlikelihood(self,nsample,query,evidence):
        closedquery=[(qt,qt) for qt in query]
        closedevidence=evidence
        openevidence=[a for (a,b) in evidence]
        querynodes,evidencenodesgiven,orderednodes,evidencenodes=self.expandNode(query,evidence)
        query=self.querytovector(orderednodes,query,evidence)
        evidence=self.querytovector(evidencenodes,query,evidence)
        for _ in range(len(query)-len(evidence)):
            evidence.append("-")
        samples=self.mldistribution(nsample,orderednodes,evidence)
        a=self.gibbsvectorcount(samples,query)
        b=self.gibbsvectorcount(samples,evidence)
        if a==0 or b==0:
            return(0)
        else:
            return(a/b)
