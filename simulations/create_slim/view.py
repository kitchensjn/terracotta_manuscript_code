import tskit
import numpy as np
import pandas as pd
import sys
sys.path.append("/Users/jameskitchens/Documents/GitHub/terracotta")
import terracotta as tct
import importlib
importlib.reload(tct)
import matplotlib.pyplot as plt


demes = pd.read_csv("dataset_glacial/demes.tsv", sep="\t")
connections = pd.read_csv("dataset_glacial/connections.tsv", sep="\t")
world_map = tct.WorldMap(demes, connections)

ts = tskit.load("/Users/jameskitchens/Desktop/slim.trees")
#random_samples = np.random.choice(ts.samples(), 50)
#ts = ts.simplify(samples=random_samples, keep_unary=True)

lineage = []
time = []
deme = []
tree = ts.first()
for sample in ts.samples():
    a = sample
    while a != -1:
        lineage.append(sample)
        time.append(ts.node(a).time)
        deme.append(int(ts.population(ts.node(a).population).metadata["name"][1:]) % len(world_map.demes))
        a = tree.parent(a)

locations = pd.DataFrame({"lineage":lineage, "time":time, "deme":deme})

previous_counts = None
for time, group in locations.groupby("time"):
    
    counts = group["deme"].value_counts()
    if (not counts.equals(previous_counts)):

        xs = []
        ys = []
        ss = []
        for i,row in world_map.demes.iterrows():
            if row["id"] in counts:
                xs.append(row["xcoord"])
                ys.append(row["ycoord"])
                ss.append(counts[row["id"]]*10)
        
        plt.scatter(world_map.demes["xcoord"], world_map.demes["ycoord"], s=1000, marker="H")
        plt.scatter(xs, ys, s=ss, color="orange", zorder=20, label="Lineages")
        
        plt.title(time)
        plt.axis("off")
        plt.show()

        previous_counts = counts


#print(np.unique_counts(ts.tables.nodes.population))