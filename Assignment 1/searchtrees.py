#Graph Filename
graphFile='graph.csv'
import sys
from collections import defaultdict

def loadGraph(graphFile):
    graph=defaultdict(list)
    file=open(graphFile,'r')
    for line in file.readlines():
        node=line.rstrip().split(',')
        graph[node[0]].append({node[1]:int(node[2])})
        graph[node[1]].append({node[0]:int(node[2])})
    return graph

def queuePop(queue,popIndex,bd):
    if bd:
        return queue.pop(popIndex)
    else:
        return queue.pop()

def goalCheck(visited,destination):
    visitedList=[k for (k,v) in visited]
    if destination in visitedList:
        return True
    else:
        return False

def resolveDistance(graph,visited):
    path=[]
    visitedDict=defaultdict(type(visited[0][0]))
    for (k,v) in visited:
        visitedDict[k]=v
    currNode=visited[len(visited)-1][0]
    while currNode is not visited[0][0]:
        path.append(currNode)
        currNode=visitedDict[currNode]
    path.append(visited[0][0])
    path.reverse()
    distance=0
    for index in range(0,len(path)-1):
        for element in graph[path[index]]:
            if list(element.keys())[0]==path[index+1]:
                distance+=element[path[index+1]]
    print(path,distance)


def bfs(source,destination,graph,bd):
    possibleLocs=graph.keys()
    if source not in possibleLocs:
        print("Source Not Found in Graph")
    elif destination not in possibleLocs:
        print("Destination Not Found in Graph")
    else:
        queue=[(source,0)]
        visited=[]
        while queue:
            currNode=queuePop(queue,0,bd)
            visited.append(currNode)
            checkGoal=goalCheck(visited,destination)
            if not checkGoal:
                for element in graph[currNode[0]]:
                    for node in element:
                        if ((node not in [k for (k,v) in visited]) and (node not in [k for (k,v) in queue])):
                            queue.append((node,currNode[0]))
            else:
                resolveDistance(graph,visited)
                return visited
        print("No Possible Goals")
        return 0

args=sys.argv[1:]
graph=loadGraph(graphFile)
bfs(args[0],args[1],graph,True)
bfs(args[0],args[1],graph,False)
