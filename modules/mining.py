#Find the parent of a given node and output a complete path.
#This function will be called and explained in the next function.
#prefixpath will be assigned as an empty list.
def completepath(node, prefixpath):
    #node is an object in class node and prefixpath is a list
    if node.parent != None:
        prefixpath.append(node.name)
        completepath(node.parent, prefixpath)

def findcompletepath(phenotype, headerTable):
    node = headerTable[phenotype][1]
    #A node class object, as the lowest level of this path, corresponding a certain phenotype.
    completeset = {}
    #a new dict to store the results: a set of all possible group of phenotypes starting from a certain phenotype.
    #(searching in one direction)
    while node != None:
        prefixpath = []
        completepath(node, prefixpath)
        #results are stored in list 'prefixpath'
        if len(prefixpath) > 1:
            completeset[frozenset(prefixpath[1:])] = node.count
        node = node.nodeLink
        #continue to the next node(but with the same phenotype)
    return completeset
    #a dict containing key : value pair -- frozenset of prefix of a phenotype :  count of this prefix.

def mineFPtree(FPTREE, headerTable, minSup, freqItemList = [], prefix = set([])):
    orderedHTB = [x[0] for x in sorted(headerTable.items(), key=lambda x:x[1])]
    for each in orderedHTB:
        newfreqset = prefix.copy()
        newfreqset.add(each)
        freqItemList.append(newfreqset)
        conditionbase = findcompletepath(each, headerTable)
        conditiontree, conditionhead = createFPtree(conditionbase, minSup)
        if conditionhead != None:
            mineFPtree(conditiontree, conditionhead, minSup, freqItemList, newfreqset)
