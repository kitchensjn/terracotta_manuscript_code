import sys
sys.path.append("/Users/jameskitchens/Documents/GitHub/terracotta")
import terracotta as tct
import importlib
importlib.reload(tct)
import pandas as pd



demes = pd.read_csv("dataset/demes.tsv", sep="\t")
connections = pd.read_csv("dataset/connections.tsv", sep="\t")
samples = pd.read_csv("dataset/samples.tsv", sep="\t")
world_map = tct.WorldMap(demes, connections, samples)

world_map.convert_to_demes_yml("demes.yaml", pop_size=10, parameters=[0.001])