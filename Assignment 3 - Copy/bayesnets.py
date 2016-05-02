from collections import defaultdict
from BayesNode import BayesNode
import copy
import pprint

pp = pprint.PrettyPrinter(indent=4)

class BayesNets:
    def __init__(self):
        """
        Class Member Declarations
        """
        self.adjacencyMetadata=defaultdict(str)
        self.adjacency=defaultdict(str)
        self.topology={'B':1,'E':2,'A':3,'J':4,'M':5}

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
        self.chainrule([],[],orderednodes,nodes)
        self.chainrule([],[],evidencenodes,evidencenodesgiven)


    def chainrule(self,query,evidence,topology,nodes):
        topstack=copy.deepcopy(topology)
        topstack.reverse()
        copystack=[]
        for i in range(len(topstack)):
            copystack.append(topstack[i:])
        trimmed=[]
        for chain in copystack:
            loctrimmed=[]
            if chain[0] in nodes:
                loctrimmed.append((chain[0],'C'))
            else:
                loctrimmed.append((chain[0],'E'))
            parents=self.adjacency[loctrimmed[0][0]].parents
            if not(parents):
                parents=[loctrimmed[0][0]]
            if len(chain)>0:
                for i in range(1,len(chain)):
                    if chain[i] in parents:
                        if chain[i] in nodes:
                            loctrimmed.append((chain[i],'C'))
                        else:
                            loctrimmed.append((chain[i],'E'))
                trimmed.append(loctrimmed)
        print(trimmed)

test=BayesNets()
test.readAdjacencyGraph("adjacencylist.txt")
test.buildAdjacencyMetdatata("cpt.txt")
test.expandNode(['B'],[('A','T')])
