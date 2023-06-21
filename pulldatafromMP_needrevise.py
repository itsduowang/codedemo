from matminer.data_retrieval.retrieve_MP import MPDataRetrieval
mpdr = MPDataRetrieval(api_key="TUvMnIDcCmqNEiG0")
from pymatgen.ext.matproj import MPRester
mpr = MPRester('TUvMnIDcCmqNEiG0')
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram
from collections import defaultdict

### testing the pourbaix diagram plot
# from pymatgen.analysis.pourbaix_diagram import PourbaixPlotter
# plotter = PourbaixPlotter(pbx)
# plotter.get_pourbaix_plot().show()

# nonmetal = ["H", "He", "B", "C", "N", "F", "Ne", "Si", "P",
#             "S", "Cl", "Ar", "As", "Se", "Br", "Kr", "Te", "I", "Xe",
#             "At", "Rn"]
# nonmetal_wO = ["H", "He", "B", "C", "N", "O", "F", "Ne", "Si", "P",
#             "S", "Cl", "Ar", "As", "Se", "Br", "Kr", "Te", "I", "Xe",
#             "At", "Rn"]
# df_MO = mpdr.get_dataframe(criteria={"nelements": 2,"elements": "$nin": nonmetal, "elements": {"$all": ["O"]}},
#                            properties=['material_id', 'pretty_formula', 'e_above_hull'])
# print("There are {} entries of metal oxides".format(df_MO['pretty_formula'].count()))

# df_BM = mpdr.get_dataframe(criteria={"nelements": 2,"elements": {"$nin": nonmetal_wO}},
#                            properties=['material_id', 'pretty_formula', 'e_above_hull'])
# print("There are {} entries of bimetallic compounds".format(df_BM['pretty_formula'].count()))

# df_TM = mpdr.get_dataframe(criteria={"nelements": 3,"elements": {"$nin": nonmetal_wO}},
#                             properties=['material_id', 'pretty_formula', 'e_above_hull'])
# print("There are {} entries of trimetallic compounds".format(df_TM['pretty_formula'].count()))

# df_MO.to_json('df_MO.json')
# df_BM.to_json('df_BM.json')
# df_TM.to_json('df_TM.json')

from pymatgen.core import Composition
def screen_metallic(mat_list, oxides=True):
    new_dict = defaultdict(dict)
    nonmetal = ["H", "He", "B", "C", "N", "F", "Ne", "Si", "P",
                "S", "Cl", "Ar", "As", "Se", "Br", "Kr", "Te", "I", "Xe",
                "At", "Rn", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
                "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Ac", "Th", "Pa", "U", "Np", "Pu"]
    nonmetal_wO = ["H", "He", "B", "C", "N", "O", "F", "Ne", "Si", "P",
                   "S", "Cl", "Ar", "As", "Se", "Br", "Kr", "Te", "I", "Xe",
                   "At", "Rn", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
                   "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Ac", "Th", "Pa", "U", "Np", "Pu"]
    screen_list = nonmetal if oxides else nonmetal_wO
    for i in range(len(mat_list)):
        adding_tag = True
        temp_key = list(mat_list.keys())[i]
        temp_value = list(mat_list.values())[i]
        MP_composition = Composition(temp_value)
        for ele in (MP_composition.elements):
            if str(ele) in screen_list:
                adding_tag = False
        if (adding_tag):
            new_dict[temp_key] = temp_value
    return new_dict

import json
f_MO = open('df_MO.json', 'r')
MO = json.load(f_MO)
MO = MO['pretty_formula']
MO_list = screen_metallic(MO)
f_BM = open('df_BM.json', 'r')
BM = json.load(f_BM)
BM = BM['pretty_formula']
BM_list = screen_metallic(BM, oxides=False)
f_TM = open('df_TM.json', 'r')
TM = json.load(f_TM)
TM = TM['pretty_formula']
TM_list = screen_metallic(TM, oxides=False)

print (len(MO),len(BM),len(TM))
print (len(MO_list),len(BM_list),len(TM_list))

def get_phxeng(mat_dict, ph_max=7, pot_mmin=-1.0, ifMO=False):
    selected_dict = defaultdict(dict)
    for d in range(len(mat_dict)):
        if d % 51 == 0:
            mpr = MPRester('TUvMnIDcCmqNEiG0')
            from time import sleep
            sleep(10)
#         if d == 20: break
        temp_key = list(mat_dict.keys())[d]
        temp_value = list(mat_dict.values())[d]
        print (d, temp_value)
        temp_composition = Composition(temp_value)
        mpcomound = []
        for i in temp_composition.elements:
            mpcomound.append(str(i))
        if (ifMO):
            mpcomound.remove('O')
            entries = mpr.get_pourbaix_entries(mpcomound)
            print (mpcomound)
            pbx = PourbaixDiagram(entries)
        else:
            entries = mpr.get_pourbaix_entries(mpcomound)
            print (len (entries))
            composition_dict = temp_composition.fractional_composition.as_dict()
            pbx = PourbaixDiagram(
                entries,
                comp_dict=composition_dict,
                filter_solids=True,
            )
#             from pymatgen.analysis.pourbaix_diagram import PourbaixPlotter
#             plotter = PourbaixPlotter(pbx)
#             plotter.get_pourbaix_plot().show()
        entry = [e for e in entries if e.entry_id == temp_key][0]
        small_phx = 100
        for pH_val in range (0, ph_max+1):
            temp_phx = pbx.get_decomposition_energy(entry, pH=pH_val, V=pot_mmin)
            if temp_phx <= small_phx:
                small_phx = temp_phx
                temp_pH = pH_val
#           print (small_phx, temp_pH)
        selected_dict[temp_key]["composition"] = temp_value
        selected_dict[temp_key]["energy"] = small_phx
        print ("done")
    return selected_dict


def split_dict(d, n):
    keys = list(d.keys())
    for i in range(0, len(keys), n):
        yield {k: d[k] for k in keys[i: i + n]}
# for item in split_dict({i: i for i in range(10)}, 3):
#     print(item)


temp_store_dict = defaultdict(dict)

breakpoint = 13
i = -1
for phx_chunk in split_dict(MO_list, 100):
    i += 1
    print (i)
    if i < breakpoint: continue
    temp_store_dict[i] = get_phxeng(phx_chunk, ifMO=True)

print (i)

for key in temp_store_dict.keys():
    print (key)


f_prev_stored_MO = open('temp_store_MO.json', 'r')
prev_stored_MO = defaultdict(dict)
prev_stored_MO = json.load(f_prev_stored_MO)
prev_stored_MO

for key in temp_store_dict.keys():
    prev_stored_MO[key] = dict(temp_store_dict[key])

with open('temp_store_MO.json', 'w') as fp:
    json.dump(prev_stored_MO, fp)

len(prev_stored_MO)


with open('store_MO.json', 'w') as fp:
    json.dump(prev_stored_MO, fp)