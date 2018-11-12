import dill

filename = 'FPTREE.txt'
txt = open(filename)
data = txt.readlines()
data = data[1:]
# Open and load the file. 'data' is a list.
disease_origin = []
phenotype_origin = []

for x in data:
    disease_i = x.split('\t')[0].strip()
    phenotype_i = x.split('\t')[1].strip()
    disease_origin.append(int(disease_i))
    phenotype_origin.append(int(phenotype_i))

disease_unique = []
phenotype_unique = []
for x in disease_origin:
    if x not in disease_unique:
        disease_unique.append(x)
for y in phenotype_origin:
    if y not in phenotype_unique:
        phenotype_unique.append(y)

convert = {}
length = len(disease_origin)
for i in range(length):
    if disease_origin[i] not in convert:
        convert[disease_origin[i]] = [phenotype_origin[i]]
    else:
        convert[disease_origin[i]].append(phenotype_origin[i])
pheno_set = list(convert.values())

class node:
    def __init__(self, name, count, parent):
        # 5 attributes.
        self.name = name
        self.count = count
        self.parent = parent
        self.nodeLink = None
        self.children = {}
        # The nodes next level to self, key : value is name : node.

    def increase(self, count):
        self.count += count
        # counter of the node

    def disp(self, ind=1):
        print('  ' * ind, self.name, ' ', self.count)
        for child in self.children.values():
            child.disp(ind + 1)
        # for showing the tree shape.

# This function turn the list of lists into a dict counting each group of phenotypes.
def createinit(pheno_set):
    initdic: Dict[FrozenSet[Any], int] = {}
    for event in pheno_set:
        key = frozenset(event)
        if key in initdic:
            initdic[key] += 1
        else:
            initdic[key] = 1
    return initdic

def updateHeader(nodeToTest, targetNode):
    while nodeToTest.nodeLink != None:
        nodeToTest = nodeToTest.nodeLink
    nodeToTest.nodeLink = targetNode

# This function will be needed in the next function.
def updateFPtree(event, inTree, headerTable, count):
    # event is a list of phenotypes of one disease,
    # inTree is a destination node, in class node.
    # headerTable is just headerTable.
    # count is the occuring number of this event. (see createinit())
    if event[0] in inTree.children:
        # check if the first item of event is already a child of inTree.
        inTree.children[event[0]].increase(count)
        # chindren is also in class node.
    else:
        # create new branch.
        inTree.children[event[0]] = node(event[0], count, inTree)
        if headerTable[event[0]][1] == None:
            headerTable[event[0]][1] = inTree.children[event[0]]
        else:
            updateHeader(headerTable[event[0]][1], inTree.children[event[0]])
    if len(event) > 1:
        updateFPtree(event[1:], inTree.children[event[0]], headerTable, count)
# The main executer of creating fptree.
def createFPtree(initdic, minSup):
    headerTable = {}
    for event in initdic:
        for phenotype in event:
            headerTable[phenotype] = headerTable.get(phenotype, 0) + initdic[event]
            # Counting for every single phenotype.
    for k in list(headerTable.keys()):
        if headerTable[k] < minSup:
            del (headerTable[k])
            # Delete the phenotypes with occuring number < minSup.
    freq_pheno_set = set(headerTable.keys())
    # Make a set of the 'freq' phenotypes.
    if len(freq_pheno_set) == 0:
        # minSup is too large.
        return None, None
    for k in headerTable:
        headerTable[k] = [headerTable[k], None]
        # Change the format: element [count, node].

    rootTree = node('Null Set', 1, None)
    # Make a new ROOT node, with name 'Null Set', count '1', and no parent.
    for event, count in initdic.items():
        # initdic：[element(group of phenotypes), count]
        select_phenotype = {}
        for phenotype in event:
            if phenotype in freq_pheno_set:
                select_phenotype[phenotype] = headerTable[phenotype][0]
                # Select the phenotypes with support > minSup.
                # In select_phenotype, phenotype : count
        if len(select_phenotype) > 0:
            # Order the phenotypes by counting number.
            ordered_list = sorted(select_phenotype.items(), key=lambda x: x[1], reverse=True)
            ordered_phenotype = []
            for each in ordered_list:
                ordered_phenotype.append(each[0])
            # Return a list of phenotype, selected and ordered.
            updateFPtree(ordered_phenotype, rootTree, headerTable, count)
    return rootTree, headerTable

initdic = createinit(pheno_set)

filehandle = open('hpo1009.txt')
fileload = filehandle.read()
txtsplit = fileload.split('\n\n')
# txtsplit is a list of chunks.
terms = txtsplit[1:]

term_id_list = []
term_id_list_inuse = []
term_name_list = []
name_dic = {}
parent_offspring = set()

import re

for term in terms:
    id = int(re.search('id: HP:(\d+)', term).group(1))
    name = re.search('\nname: (.+)\n', term).group(1)
    parents = re.findall('is_a: HP:(\d+)', term)
    if len(parents) != 0:
        term_id_list_inuse.append(id)
        for i in range(len(parents)):
            parents[i] = int(parents[i])
    term_id_list.append(id)
    term_name_list.append(name)
    if len(parents) != 0:
        for i in range(len(parents)):
            parent_offspring.add((parents[i], id))

for i in range(len(term_id_list)):
    name_dic[term_id_list[i]] = term_name_list[i]


import networkx as nx

from networkx import all_neighbors, ancestors, descendants
list_PO = list(parent_offspring)
G = nx.DiGraph()
G.add_edges_from(list_PO)
# convert is a dict ,{disease_id : list of specific phenotypes}
# print(convert)
# Inflate the dict to calculate information content.
# print(convert.items())
def inflate(ls):
    for each in ls:
        ls.extend(list(ancestors(G, each)))
    return set(ls)

convert_inflate = convert
for each in list(convert_inflate.keys()):
    convert_inflate[each] = inflate(convert_inflate[each])

AnnoNum = {}
# Dict, phenoID : number of annotation
for each in convert_inflate.values():
    each = list(each)
    for every in each:
        if every not in AnnoNum:
            AnnoNum[every] = 1
        if every in AnnoNum:
            AnnoNum[every] += 1

dillfile = 'base.pkl'
dill.dump_session(dillfile)
