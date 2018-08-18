#This is for creating the fptree structure.
filename = 'FPTREE.txt'
txt = open(filename)
data = txt.readlines()
#Open and load the file. 'data' is a list.

disease = []
phenotype = []
for x in data:
    disease_i = x.split('\t')[0].strip()
    phenotype_i = x.split('\t')[1].strip()
    disease.append(disease_i)
    phenotype.append(phenotype_i)
#Delete the header of the table.
disease = disease[1:]
phenotype = phenotype[1:]
#disease and phenotype are lists now.

convert = {}
length = len(disease)
#length of disease and phenotype are equal.
for i in range(length):
    if disease[i] not in convert:
        convert[disease[i]] = [phenotype[i]]
    else:
        convert[disease[i]].append(phenotype[i])
#convert is a dict, the key : value is disease_id : phenotype_id(in list).
pheno_set = list(convert.values())
#A list of lists, each element list is a group of phenotype occuring together.
#pheno_set is ready to go in createinit()

#Define a class called node.
class node:
    def __init__(self, name, count, parent):
        #5 attributes.
        self.name = name
        self.count = count
        self.parent = parent
        self.nodeLink = None
        self.children = {}
        #The nodes next level to self, key : value is name : node.

    def increase(self, count):
        self.count += count
        #counter of the node

    def disp(self, ind=1):
        print( '  '*ind, self.name, ' ', self.count)
        for child in self.children.values():
            child.disp(ind+1)
        #for showing the tree shape.

#This function turn the list of lists into a dict counting each group of phenotypes.
def createinit(pheno_set):
    initdic = {}
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



#This function will be needed in the next function.
def updateFPtree(event, inTree, headerTable, count):
    #event is a list of phenotypes of one disease,
    #inTree is a destination node, in class node.
    #headerTable is just headerTable.
    #count is the occuring number of this event. (see createinit())
    if event[0] in inTree.children:
        #check if the first item of event is already a child of inTree.
        inTree.children[event[0]].increase(count)
        #chindren is also in class node.
    else:
        #create new branch.
        inTree.children[event[0]] = node(event[0], count, inTree)
        if headerTable[event[0]][1] == None:
            headerTable[event[0]][1] = inTree.children[event[0]]
        else:
            updateHeader(headerTable[event[0]][1], inTree.children[event[0]])
    if len(event) > 1:
        updateFPtree(event[1:], inTree.children[event[0]], headerTable, count)

#The main executer of creating fptree.
def createFPtree(initdic, minSup=1):
    headerTable = {}
    for event in initdic:
        for phenotype in event:
            headerTable[phenotype] = headerTable.get(phenotype, 0) + initdic[event]
            #Counting for every single phenotype.
    for k in list(headerTable.keys()):
        if headerTable[k] < minSup:
            del(headerTable[k])
            #Delete the phenotypes with occuring number < minSup.
    freq_pheno_set = set(headerTable.keys())
    #Make a set of the 'freq' phenotypes.
    if len(freq_pheno_set) == 0:
        #minSup is too large.
        return None, None
    for k in headerTable:
        headerTable[k] = [headerTable[k], None]
        #Change the format: element [count, node].

    rootTree = node('Null Set', 1, None)
    #Make a new ROOT node, with name 'Null Set', count '1', and no parent.
    for event, count in initdic.items():
        # initdic：[element(group of phenotypes), count]
        select_phenotype = {}
        for phenotype in event:
            if phenotype in freq_pheno_set:
                select_phenotype[phenotype] = headerTable[phenotype][0]
                #Select the phenotypes with support > minSup.
                #In select_phenotype, phenotype : count
        if len(select_phenotype) > 0:
            #Order the phenotypes by counting number.
            ordered_list = sorted(select_phenotype.items(), key=lambda x:x[1], reverse=True)
            ordered_phenotype = []
            for each in ordered_list:
                ordered_phenotype.append(each[0])
            #Return a list of phenotype, selected and ordered.
            updateFPtree(ordered_phenotype, rootTree, headerTable, count)
    return rootTree, headerTable
#Reference:
#https://blog.csdn.net/songbinxu/article/details/80411388
