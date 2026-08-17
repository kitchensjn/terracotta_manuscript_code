import sys
sys.path.append("/Users/jameskitchens/Documents/GitHub/terracotta")
import terracotta as tct
import importlib
importlib.reload(tct)
sys.path.append("/Users/jameskitchens/Documents/GitHub/terracotta_manuscript_code/simulations")
import sim_utils
import pandas as pd


sim_utils.create_hexagonal_tri_grid(2, output_directory="dataset_glacial_small")

demes = pd.read_csv("dataset_glacial_small/demes.tsv", sep="\t")
demes["suitability"] = demes["suitability"].astype(str)
demes.loc[demes["ycoord"] > 0, "suitability"] = "0:1,500:0"

demes.to_csv("dataset_glacial_small/demes.tsv", index=False, sep="\t")