import sys
import BayesNets
args=sys.argv[1:]

bayesnet=BayesNets.BayesNets()
bayesnet.readAdjacencyGraph("adjacencylist.txt")
bayesnet.buildAdjacencyMetdatata("cpt.txt")

def getinput():
    line_holders=input().split(" ")
    query_variables=[]
    evidence_variables=[]
    line=0
    for line in range(int(line_holders[0])):
        evidence_variables.append(tuple([x.upper() for x in input().split(" ")]))
    for line in range(int(line_holders[1])):
        query_variables.append(input())
    return(query_variables,evidence_variables)

def evaluate(q,e):
    o=10
    if args[0]=='e':
        a=[]
        for item in q:
            a=[]
            for _ in range(o):
                a.append(bayesnet.enumeration([item],e))
            print((item,sum(a)/len(a)))
            #print(bayesnet.enumeration(['B'],[('M','T'),('J','T')]))
    elif args[0]=='p':
        a=[]
        for item in q:
            a=[]
            for _ in range(o):
                a.append(bayesnet.priorsampling(float(args[1]),item,e))
            print((item,sum(a)/len(a)))
            #print(bayesnet.priorsampling(float(args[1]),['J'],[('A','T'),('B','F')]))
    elif args[0]=='r':
        a=[]
        for item in q:
            a=[]
            for _ in range(o):
                a.append(bayesnet.rejectionsampling(float(args[1]),item,e))
            print((item,sum(a)/len(a)))
        #print(bayesnet.priorsampling(float(args[1]),['J'],[('A','T'),('B','F')]))
    elif args[0]=='l':
        for item in q:
            a=[]
            for _ in range(o):
                a.append(bayesnet.maxlikelihood(float(args[1]),item,e))
            print((item,sum(a)/len(a)))
        #print(bayesnet.priorsampling(float(args[1]),['J'],[('A','T'),('B','F')]))
    else:
        print("Invalid Parameters")

q,e=getinput()
for k in [10,50,100,200,500,1000,10000]:
    for j in ['p','r','l','e']:
        args[0]=j
        args[1]=k
        print("-------------------------------------------")
        print(args[0],args[1])
        print("-------------------------------------------")
        evaluate(q,e)
