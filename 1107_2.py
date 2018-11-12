import datetime
starttime = datetime.datetime.now()
import networkx as nx

from networkx import all_neighbors, ancestors, descendants
import dill
dillfile = 'base.pkl'
dill.load_session(dillfile)

AnnoNum = {}
# Dict, phenoID : number of annotation
for each in convert_inflate.values():
    each = list(each)
    for every in each:
        if every not in AnnoNum:
            AnnoNum[every] = 1
        if every in AnnoNum:
            AnnoNum[every] += 1

print(AnnoNum)
import math
def GetIC(term):
    return (math.log(AnnoNum[1]/AnnoNum[term]))

def GetMICA_IC(term1, term2):
    ancestor_set1 = ancestors(G, term1)
    ancestor_set2 = ancestors(G, term2)
    common_set = ancestor_set1 & ancestor_set2
    sorted_list = sorted(list(common_set), key = lambda x: GetIC(x), reverse=True)
    return GetIC(sorted_list[0])



print(GetIC(1))
print(GetIC(10))
print(GetMICA_IC(44,8))
endtime = datetime.datetime.now()
print (endtime - starttime)