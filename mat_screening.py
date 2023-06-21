import json
import pprint
from collections import defaultdict
from pymatgen.core.structure import *
from matminer.data_retrieval.retrieve_MP import MPDataRetrieval

f_MO = open('store_MO.json', 'r')
MO_0 = json.load(f_MO)
f_BM = open('store_BM.json', 'r')
BM_0 = json.load(f_BM)
f_MO = open('store_TM.json', 'r')
TM_0 = json.load(f_MO)


def full_list(list_0):
    list_all = defaultdict(dict)
    count = 0
    for i in list_0.values():
        #     print (len(i))
        for j in list(i.keys()):
            #           print (i[j]['composition'])
            list_all[count]['mat_id'] = j
            list_all[count]['composition'] = i[j]['composition']
            list_all[count]['pbxenergy'] = i[j]['energy']
            count += 1
    return list_all


### filter 1: PBX
### filter 2: ehull
### filter 3: cost
def pbxfilter(list_in, energy=0.2):
    list_out = defaultdict(dict)
    count = 0
    for i in range(len(list_in)):
        if list_in[i]['pbxenergy'] < energy:
            list_out[count] = list_in[i]
            count += 1
    return list_out


from pymatgen.analysis.cost import CostAnalyzer, CostDBCSV

costdb_path = "/Users/duowang/Documents/GitHub/nitrate/datasets/costdb_elements_2021.csv"
costdb = CostDBCSV(costdb_path)
costanalyzer = CostAnalyzer(costdb)


def ehullfilter(list_in, energy_barrier=0.1):
    mpdr = MPDataRetrieval(api_key="TUvMnIDcCmqNEiG0")
    list_out = defaultdict(dict)
    count = 0
    mpidlist = [];
    for i in range(len(list_in)):
        mpidlist.append(list_in[i].get('mat_id'))
    #    print (mpidlist)
    temp_data = mpdr.get_dataframe({"material_id": {"$in": mpidlist}},
                                   properties=["material_id", 'pretty_formula', "e_above_hull",
                                               'formation_energy_per_atom'])
    #     print (temp_data)
    temp_data.to_json('temp_data.json')
    temp_file = open('temp_data.json', 'r')
    temp_data = json.load(temp_file)
    #     print (list(temp_data.values())[2])
    dict_idtocom = list(temp_data.values())[0]
    dict_idtoene = list(temp_data.values())[1]
    dict_idtoforene = list(temp_data.values())[2]
    for i in range(len(dict_idtocom)):
        temp_matid = list(dict_idtoene.keys())[i]
        temp_comp = list(dict_idtocom.values())[i]
        temp_energy = list(dict_idtoene.values())[i]
        temp_forene = list(dict_idtoforene.values())[i]
        if temp_energy <= energy_barrier:
            list_out[count]['mat_id'] = temp_matid
            list_out[count]['composition'] = temp_comp
            list_out[count]['ehull'] = temp_energy
            list_out[count]['formation_ene'] = temp_forene
            count += 1
    return (list_out)


def costfilter(list_in, cost=3000):
    list_out = defaultdict(dict)
    count = 0
    for i in range(len(list_in)):
        comp = Composition(list_in[i]['composition']).reduced_formula
        temp_cost = costanalyzer.get_cost_per_kg(comp)
        if temp_cost < cost:
            list_out[count] = list_in[i]
            list_out[count]['unitprice'] = temp_cost
            count += 1
    return list_out


def rmDupDict(list_in):
    list_out = defaultdict(dict)
    seen = []
    count = 0
    for i in range(len(list_in)):
        comp = Composition(list_in[i]['composition']).reduced_formula
        if comp not in seen:
            list_out[count] = list_in[i]
            seen.append(comp)
            count += 1
        if comp in seen:
            if list_in[i]['formation_ene'] < list_out[seen.index(comp)]['formation_ene']:
                list_out[seen.index(comp)] = list_in[i]
    return list_out


def get_compounds(dict_in):
    compound_list = [];
    for i in range(len(dict_in)):
        compound_list.append(dict_in[i].get('composition'))
    return compound_list


def rmDup(list_in):
    list_out = []
    for i in list_in:
        if i not in list_out:
            list_out.append(i)
    return list_out


def includeElement(dict_in, ele_list=['Ni']):
    dict_out = defaultdict(dict)
    count = 0
    for i in range(len(dict_in)):
        comp = Composition(dict_in[i]['composition']).reduced_formula
        tag = 0
        for j in ele_list:
            if j in comp:
                tag = 1
        if tag == 1:
            dict_out[count] = dict_in[i]
            count += 1
    return dict_out


def notincludeElement(dict_in, ele_list=['As', 'Be', 'Cd', 'Cr', 'Pb', 'Hg']):
    dict_out = defaultdict(dict)
    count = 0
    for i in range(len(dict_in)):
        comp = Composition(dict_in[i]['composition']).reduced_formula
        tag = 0
        for j in ele_list:
            if j in comp:
                tag = 1
        if tag == 0:
            dict_out[count] = dict_in[i]
            count += 1
    return dict_out


def get_pairs(list_in):
    pair_list = []
    for item in list_in:
        comp = Composition(item).as_dict()
        p = tuple([c[0] for c in sorted(list(comp.items()), reverse=True, key=lambda c: c[1])])
        if p not in pair_list:
            pair_list.append(p)
    return sorted(pair_list, key=lambda x: x[0])


def rmDupPair(plist_in):
    plist_out, seen = [], set()
    for item in plist_in:
        t1 = tuple(item)
        if t1 not in seen and tuple(reversed(item)) not in seen:
            seen.add(t1)
            plist_out.append(item)
    return sorted(plist_out, key=lambda x: x[0])


def tuple3to2(tlist3):
    t1 = [x[:2] for x in tlist3]
    t2 = [x[1:3] for x in tlist3]
    t3 = [x[0:3:2] for x in tlist3]
    t_all = t1 + t2 + t3
    return t_all


def printTable(plist_in):
    out_list = []
    for i in plist_in:
        if i[0] not in out_list:
            out_list.append(i[0])
            print(i[0], end=": ")
            seen = []
            for j in plist_in:
                if (j[0] == i[0]):
                    if j[1] not in seen:
                        seen.append(j[1])
                        print(j[1], end=" ")
            print('\n')


def getMOElements(plist_in):
    colA = []
    for i in MO_pairs:
        if i[0] != 'O' and i[0] not in colA:
            colA.append(i[0])
        elif i[0] not in colA:
            colA.append(i[1])
    return colA

MO_all = full_list(MO_0)
BM_all = full_list(BM_0)
TM_all = full_list(TM_0)
print (len(MO_all),len(BM_all),len(TM_all))
MO_filter1 = pbxfilter(MO_all)
BM_filter1 = pbxfilter(BM_all)
TM_filter1 = pbxfilter(TM_all)
print (len(MO_filter1),len(BM_filter1),len(TM_filter1))
MO_filter2 = ehullfilter(MO_filter1)
BM_filter2 = ehullfilter(BM_filter1)
TM_filter2 = ehullfilter(TM_filter1)
print (len(MO_filter2),len(BM_filter2),len(TM_filter2))
MO_filter3 = costfilter(MO_filter2, cost=3000)
BM_filter3 = costfilter(BM_filter2, cost=3000)
TM_filter3 = costfilter(TM_filter2, cost=3000)
print (len(MO_filter3),len(BM_filter3),len(TM_filter3))
# with open('filtered_MO500.json', 'w') as op_MO:
#     json.dump(MO_filter3, op_MO)
# with open('filtered_BM500.json', 'w') as op_BM:
#     json.dump(BM_filter3, op_BM)
# with open('filtered_TM500.json', 'w') as op_TM:
#     json.dump(TM_filter3, op_TM)


MO_noDup = rmDupDict(MO_filter3)
print (len(MO_noDup))

BM_noDup = rmDupDict(BM_filter3)
print (len(BM_noDup))

TM_noDup = rmDupDict(TM_filter3)
print (len(TM_noDup))

# with open('filtered_MO.json', 'w') as op_MO:
#     json.dump(MO_noDup, op_MO)
# with open('filtered_BM.json', 'w') as op_BM:
#     json.dump(BM_noDup, op_BM)
# with open('filtered_TM.json', 'w') as op_TM:
#     json.dump(TM_noDup, op_TM)


ele_remove_list=['As', 'Be', 'Cd', 'Cr', 'Pb', 'Hg']
BM_selected = notincludeElement(BM_noDup, ele_remove_list)
# BM_selected = rmDupDict(BM_selected)
print (len(BM_selected))

ele_keep_list = ['Ag', 'Cu', 'Fe', 'Ni']
BM_selected = includeElement(BM_selected, ele_keep_list)
# BM_selected = rmDupDict(BM_selected)
print (len(BM_selected))
print (BM_selected)




# MO_compounds = get_compounds(MO_filter3)
# MO_pairs = get_pairs(MO_compounds)
# eleList = getMOElements(MO_pairs)

# MO_compounds = get_compounds(MO_filter3)
# MO_comp_nodup = rmDup(MO_compounds)
# BM_compounds = get_compounds(BM_filter3)
# BM_comp_nodup = rmDup(BM_compounds)
# TM_compounds = get_compounds(TM_filter3)
# TM_comp_nodup = rmDup(TM_compounds)
# print (len(MO_comp_nodup), len(BM_comp_nodup), len(TM_comp_nodup))
# BM_pairs = get_pairs(BM_compounds)
# BM_pairs = rmDupPair(BM_pairs)
# #printTable(BM_pairs)

# TM_compounds = get_compounds(TM_filter3)
# TM_pairs3 = get_pairs(TM_compounds)
# TM_pairs = tuple3to2(TM_pairs3)
# TM_pairs = rmDupPair(TM_pairs)
# printTable(TM_pairs)

