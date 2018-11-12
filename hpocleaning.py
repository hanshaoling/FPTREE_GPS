filehandle = open('hpo1009.txt')
fileload = filehandle.read()
txtsplit = fileload.split('\n\n')
# txtsplit is a list of chunks.
terms = txtsplit[1:]

term_id_list = []
term_name_list = []
name_dic = {}
offspring_parent = set()

import re

for term in terms:
    id = int(re.search('id: HP:(\d+)', term).group(1))
    name = re.search('\nname: (.+)\n', term).group(1)
    parents = re.findall('is_a: HP:(\d+)', term)
    for i in range(len(parents)):
        parents[i] = int(parents[i])
    term_id_list.append(id)
    term_name_list.append(name)
    for i in range(len(parents)):
        offspring_parent.add((id, parents[i]))

for i in range(len(term_id_list)):
    name_dic[term_id_list[i]] = term_name_list[i]
