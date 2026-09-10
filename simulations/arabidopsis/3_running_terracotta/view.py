import sys
sys.path.append("/Users/jameskitchens/Documents/GitHub/terracotta")
import terracotta as tct
import pandas as pd
import tskit


demes = pd.read_csv("../2_creating_input_files/dataset/demes.tsv", sep="\t")
connections = pd.read_csv("../2_creating_input_files/dataset/connections.tsv", sep="\t")
samples = pd.read_csv("../2_creating_input_files/dataset/samples.tsv", sep="\t")
world_map = tct.WorldMap(demes, connections, samples)

ts = tskit.load("../2_creating_input_files/dataset/trees/0.trees")
tree = ts.simplify()

print(tree.draw_text())


exit()

print(tct.track_lineage_over_time(
    sample=0,
    times=[0,100,200,300],
    tree=tree,
    world_map=world_map,
    parameters=[0.00036968]
))