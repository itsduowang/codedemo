from matminer.data_retrieval.retrieve_MP import MPDataRetrieval
mpdr = MPDataRetrieval(api_key="TUvMnIDcCmqNEiG0")
from pymatgen.ext.matproj import MPRester
mpr = MPRester('TUvMnIDcCmqNEiG0')
from pymatgen.core import Composition

nonmetal = ["H", "He", "B", "C", "N", "F", "Ne", "Si", "P",
            "S", "Cl", "Ar", "As", "Se", "Br", "Kr", "Te", "I", "Xe",
            "At", "Rn"]
nonmetal_wO = ["H", "He", "B", "C", "N", "O", "F", "Ne", "Si", "P",
            "S", "Cl", "Ar", "As", "Se", "Br", "Kr", "Te", "I", "Xe",
            "At", "Rn"]
df_MO = mpdr.get_dataframe(criteria={"nelements": 2,"elements": {"$nin": nonmetal}, "elements": {"$all": ["O"]}},
                           properties=['material_id', 'pretty_formula', 'e_above_hull'])
print("There are {} entries of metal oxides".format(df_MO['pretty_formula'].count()))

df_BM = mpdr.get_dataframe(criteria={"nelements": 2,"elements": {"$nin": nonmetal_wO}},
                           properties=['material_id', 'pretty_formula', 'e_above_hull'])
print("There are {} entries of bimetallic compounds".format(df_BM['pretty_formula'].count()))

df_TM = mpdr.get_dataframe(criteria={"nelements": 3,"elements": {"$nin": nonmetal_wO}},
                            properties=['material_id', 'pretty_formula', 'e_above_hull'])
print("There are {} entries of trimetallic compounds".format(df_TM['pretty_formula'].count()))

df_MO.to_json('df_MO.json')
df_BM.to_json('df_BM.json')
df_TM.to_json('df_TM.json')
def screen_metallic(mat_list, oxides=True):
    new_dict = defaultdict(dict)
    nonmetal = ["H", "He", "B", "C", "N", "F", "Ne", "Si", "P",
            "S", "Cl", "Ar", "As", "Se", "Br", "Kr", "Te", "I", "Xe",
            "At", "Rn"]
    nonmetal_wO = ["H", "He", "B", "C", "N", "O", "F", "Ne", "Si", "P",
            "S", "Cl", "Ar", "As", "Se", "Br", "Kr", "Te", "I", "Xe",
            "At", "Rn"]
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


def get_phxeng(mat_dict, ph_max = 7, pot_mmin = -1.0):
    selected_dict = defaultdict(dict)
    for d in range(len(mat_dict)):
        if d % 10 == 0: mpr = MPRester('TUvMnIDcCmqNEiG0')
        temp_key = list(mat_dict.keys())[d]
        temp_value = list(mat_dict.values())[d]
        temp_composition = Composition(temp_value)
        mpcomound = []
        for i in temp_composition.elements:
            mpcomound.append(str(i))
        entries = mpr.get_pourbaix_entries(mpcomound)
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
        selected_dict[i]["material_id"] = temp_key
        selected_dict[i]["composition"] = temp_value
        selected_dict[i]["phx_energy"] = small_phx
    return selected_dict


phx_MO = get_phxeng(BM_list)
phx_BM = get_phxeng(BM_list)
phx_TM = get_phxeng(TM_list)