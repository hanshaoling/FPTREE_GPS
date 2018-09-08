# This is for creating the fptree structure.
from typing import FrozenSet, Dict, Any

filename = 'FPTREE.txt'
txt = open(filename)
data = txt.readlines()
# Open and load the file. 'data' is a list.

disease_origin = []
phenotype_origin = []
for x in data:
    disease_i = x.split('\t')[0].strip()
    phenotype_i = x.split('\t')[1].strip()
    disease_origin.append(disease_i)
    phenotype_origin.append(phenotype_i)
# Delete the header of the table.
disease_origin = disease_origin[1:]
phenotype_origin = phenotype_origin[1:]

disease_dict = {}
# pairing origin id and integer
disease_unique = []
for x in disease_origin:
    if x not in disease_unique:
        disease_unique.append(x)

i = 1
for x in disease_unique:
    if x not in disease_dict:
        disease_dict[x] = i
    i += 1

phenotype_dict = {}
phenotype_unique = []
for x in phenotype_origin:
    if x not in phenotype_unique:
        phenotype_unique.append(x)

i = 1
for x in phenotype_unique:
    if x not in phenotype_dict:
        phenotype_dict[x] = i
    i += 1

disease = []
for x in disease_origin:
    disease.append(disease_dict[x])
phenotype = []
for x in phenotype_origin:
    phenotype.append(phenotype_dict[x])
# disease and phenotype are lists now, using
# integer to store information, instead of string.

convert = {}
length = len(disease)
# length of disease and phenotype are equal.
for i in range(length):
    if disease[i] not in convert:
        convert[disease[i]] = [phenotype[i]]
    else:
        convert[disease[i]].append(phenotype[i])
# convert is a dict, the key : value is disease_id : phenotype_id(in list).
pheno_set = list(convert.values())
# A list of lists, each element list is a group of phenotype occuring together.

# Define a class called node.
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
def createFPtree(initdic, minSup = 10):
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
threshold = input('The minimal SUPPORT is set to: ')
if threshold == '':
    print('Default value: 10')
    threshold = 10
else:
    threshold = int(threshold)
FPTREE, HEADERTABLE = createFPtree(initdic, threshold)

# Find the parent of a given node and output a complete path.
# This function will be called and explained in the next function.
# prefixpath will be assigned as an empty list.
def completepath(node, prefixpath):
    # node is an object in class node and prefixpath is a list
    if node.parent != None:
        prefixpath.append(node.name)
        completepath(node.parent, prefixpath)

def findcompletepath(phenotype, headerTable):
    node = headerTable[phenotype][1]
    # A node class object, as the lowest level of this path, corresponding a certain phenotype.
    completeset = {}
    # a new dict to store the results: a set of all possible group of phenotypes starting from a certain phenotype.
    # (searching in one direction)
    while node != None:
        prefixpath = []
        completepath(node, prefixpath)
        # results are stored in list 'prefixpath'
        if len(prefixpath) > 1:
            completeset[frozenset(prefixpath[1:])] = node.count
        node = node.nodeLink
        # continue to the next node(but with the same phenotype)
    return completeset
    # a dict containing key : value pair -- frozenset of prefix of a phenotype :  count of this prefix.

# Capsule these functions.
def mineFPtree(FPTREE, headerTable, minSup = 10, freqItemList = [], prefix = set([])):
    ordered = sorted(headerTable.items(), key=lambda x: x[0])
    orderedHTB = []
    for each in ordered:
        orderedHTB.append(each[0])

    for each in orderedHTB:
        newfreqset = prefix.copy()
        newfreqset.add(each)
        freqItemList.append(newfreqset)
        conditionbase = findcompletepath(each, headerTable)
        conditiontree, conditionhead = createFPtree(conditionbase, minSup)
        if conditionhead != None:
            mineFPtree(conditiontree, conditionhead, minSup, freqItemList, newfreqset)
    return freqItemList

# For reverse rearching in dict.
def reverse_lookup(dict, value):
    for key in dict:
        if dict[key] == value:
            return key

final = mineFPtree(FPTREE, HEADERTABLE, threshold)
# Sort the list by length.
sorted_list = sorted(final, key=lambda x: len(x), reverse=True)

# Get all REAL subset of a listed-freqset. Return is a list of list. [[...], [...], [...], ...]
def getsubset(set):
    N = len(set)
    subset = []
    for i in range(2 ** N):
        combo = []
        for j in range(N):
            if ( i >> j ) % 2 == 1:
                combo.append(set[j])
        if combo != [] and len(combo) != N:
            subset.append(combo)
    return subset

# Counting for subset.
memo = {}
def count_subset(subset): # subset is a frozenset derived from list.
    ls = list(initdic.keys())
    count = 0
    if subset in memo:
        return memo[subset]
    for each in ls:
        if subset.issubset(each):
            count += 1
    memo[subset] = count
    return count

def get_association(sorted_list, minconfidence = 0.5):
    rules = {}
    for freqset in sorted_list:
        freq_listed = list(freqset)
        freq_subset = getsubset(freq_listed)
        for each in freq_subset:
            a = frozenset(each)
            count = count_subset(a)
            confidence  = count_subset(frozenset(freqset)) / count
            if confidence >= minconfidence:
                goodrule = (a, frozenset(freqset) - a)
                rules[goodrule] = confidence
    return rules

import re
filehandle = open('hp.obo')
fileload = filehandle.read()
txtsplit = fileload.split('\n\n')
# txtsplit is a list of chunks.
terms = txtsplit[2:]

# Delete the first element: description of the GO, and the second element: ROOT of the ontology.
del terms[-1]
# Delete the final element: it is a NONE.
term_id_list = []
term_name_list = []
name_dic = {}

for term in terms:
    id = int(re.search('id: HP:(\d+)', term).group(1))
    name = re.search('\nname: (.+)\n', term).group(1)
    term_id_list.append(id)
    term_name_list.append(name)

for i in range(len(term_id_list)) :
    name_dic[term_id_list[i]] = term_name_list[i]

# Convert the output from index to origin phenotype id.
sorted_list_origin = []
for i in range(len(sorted_list)):
    sorted_list_origin.append(list(sorted_list[i]))
    for j in range(len(sorted_list[i])):
        sorted_list_origin[i][j] = reverse_lookup(phenotype_dict, list(sorted_list[i])[j])
        sorted_list_origin[i][j] = name_dic[int(sorted_list_origin[i][j])]
sorted_list_origin = sorted(sorted_list_origin, key=lambda x: len(x), reverse=True)

print('The frequent sets: ')
print(*sorted_list_origin, sep = '\n')
print('*' * 20)

confidence_threshold = input('The minimal CONFIDENCE is set to: ')
if confidence_threshold == '':
    print('Default value: 0.5')
    confidence_threshold = 0.5
else:
    confidence_threshold = float(confidence_threshold)
rules = get_association(sorted_list, confidence_threshold)

def lift_filter(FROM, TO):
    n_total = len(disease_unique)
    lift = (memo[frozenset(list(FROM) + list(TO))] * n_total) / (memo[FROM] * memo[TO])
    if lift < 1:
        return 'Conflicting! LIFT = %f' %(lift)
    elif lift == 1:
        return 'Independent! LIFT = %f' %(lift)
    elif 1 < lift < 3:
        return 'Valid! LIFT = %f' %(lift)
    else:
        return 'Strong valid! LIFT = %f' %(lift)

for key, value in rules.items():
    key_list = list(key)
    lift_each = lift_filter(key_list[0], key_list[1])
    for i in range(len(key_list)):
        key_list[i] = list(key_list[i])
        for j in range(len(key_list[i])):
            key_list[i][j] = reverse_lookup(phenotype_dict, key_list[i][j])
            key_list[i][j] = name_dic[int(key_list[i][j])]
    print('Rules: ', key_list[0], 'to', key_list[1], '| Confidence: ', value, '|', lift_each)
