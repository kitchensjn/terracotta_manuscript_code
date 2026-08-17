import tskit

ts = tskit.load("dataset/trees/0.trees")
print(ts.tables.nodes[-1].time)