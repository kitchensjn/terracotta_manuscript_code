import sys
sys.path.append("/Users/jameskitchens/Documents/GitHub/terracotta")
import terracotta as tct
import importlib
importlib.reload(tct)
import pandas as pd
import numpy as np



demes = pd.read_csv("dataset_glacial_small/demes.tsv", sep="\t")
connections = pd.read_csv("dataset_glacial_small/connections.tsv", sep="\t")
world_map = tct.WorldMap(demes, connections)

start_time_in_past = 1000
default_pop_size = 15
default_migration_rate = 0.001

current_counter = [0 for i in world_map.demes["id"]]
current_deme_ids = [None for i in world_map.demes["id"]]
current_suitabilities = [0 for i in range(len(world_map.demes))]
current_start_time = start_time_in_past - 1

with open("output.txt", "w") as output:
    for epoch in world_map.epochs[::-1]:
        new_demes = []
        suitabilities = world_map.get_all_deme_suitabilities_at_time(epoch)
        output.write(str(int(start_time_in_past - current_start_time)) + " early() {\n")
        for i, deme in world_map.demes.iterrows():
            deme_size = int(suitabilities[i]*default_pop_size)  
            if current_suitabilities[i] < 1:
                if deme_size >= 1:
                    id = str(i + len(world_map.demes["id"]) * current_counter[i])
                    output.write('\tsim.addSubpop("p' + id + '", ' + str(deme_size) + ');\n')
                    new_demes.append(i)
                    current_counter[i] += 1
                    current_deme_ids[i] = id
            else:
                output.write("\tp" + str(i) + ".setSubpopulationSize(" + str(deme_size) + ");\n")
                if deme_size < 1:
                    current_deme_ids[i] = None
        if epoch != world_map.epochs[-1]:
            if len(new_demes) > 0:
                output.write("\n")
            for i in new_demes:
                origin = list(world_map.connections["deme_0"][np.where(world_map.connections["deme_1"] == world_map.demes["id"][i])[0]]) + list(world_map.connections["deme_1"][np.where(world_map.connections["deme_0"] == world_map.demes["id"][i])[0]])
                labels = ["p"+current_deme_ids[o] for o in origin if (o not in new_demes) and (current_deme_ids[o] is not None)]
                ratios = [str(1/len(labels)) for i in range(len(labels))]
                output.write("\tp" + current_deme_ids[i] + ".setMigrationRates(c(" + ",".join(labels) + "), c(" + ",".join(ratios) + "));\n")
        for i, connection in world_map.connections.iterrows():
            d0 = np.where(connection["deme_0"] == world_map.demes["id"])[0][0]
            d1 = np.where(connection["deme_1"] == world_map.demes["id"])[0][0]
            if epoch == world_map.epochs[-1]:
                d0_label = current_deme_ids[d0]
                d1_label = current_deme_ids[d1]
                if (d0_label is not None) and (d1_label is not None):
                    output.write("\tp" + d0_label + ".setMigrationRates(p" + d1_label + ", " + str(default_migration_rate) + ");\n")
                    output.write("\tp" + d1_label + ".setMigrationRates(p" + d0_label + ", " + str(default_migration_rate) + ");\n")
            elif (d0 not in new_demes) and (d1 not in new_demes):
                d0_label = current_deme_ids[d0]
                d1_label = current_deme_ids[d1]
                if (d0_label is not None) and (d1_label is not None):
                    output.write("\tp" + d0_label + ".setMigrationRates(p" + d1_label + ", " + str(default_migration_rate) + ");\n")
                    output.write("\tp" + d1_label + ".setMigrationRates(p" + d0_label + ", " + str(default_migration_rate) + ");\n")
        output.write("}\n\n")
        if epoch != world_map.epochs[-1]:
            if len(new_demes) > 0:
                output.write(str(int(start_time_in_past - current_start_time + 1)) + " early() {\n")
                for i, connection in world_map.connections.iterrows():
                    d0_label = current_deme_ids[np.where(connection["deme_0"] == world_map.demes["id"])[0][0]]
                    d1_label = current_deme_ids[np.where(connection["deme_1"] == world_map.demes["id"])[0][0]]
                    if (d0_label is not None) and (d1_label is not None):
                        output.write("\tp" + d0_label + ".setMigrationRates(p" + d1_label + ", " + str(default_migration_rate) + ");\n")
                        output.write("\tp" + d1_label + ".setMigrationRates(p" + d0_label + ", " + str(default_migration_rate) + ");\n")
                output.write("}\n\n")
        current_suitabilities = suitabilities.copy()
        current_start_time = epoch


exit()

with open("output.txt", "w") as output:
    for epoch in world_map.epochs[::-1]:
        suitabilities = world_map.get_all_deme_suitabilities_at_time(epoch)
        output.write(str(int(start_time_in_past - current_start_time)) + " early() {\n")
        for i, deme in world_map.demes.iterrows():
            deme_size = int(suitabilities[i]*default_pop_size)  
            if current_suitabilities[i] < 1:
                if deme_size >= 1:
                    id = str(i + len(world_map.demes["id"]) * current_counter[i])
                    output.write('\tsim.addSubpop("p' + id + '", ' + str(deme_size) + ');\n')
                    current_counter[i] += 1
                    current_deme_ids[i] = id
            else:
                output.write("\tp" + str(i) + ".setSubpopulationSize(" + str(deme_size) + ");\n")
                if deme_size < 1:
                    current_deme_ids[i] = None
        output.write("\n")
        for i, connection in world_map.connections.iterrows():
            d0_label = current_deme_ids[np.where(connection["deme_0"] == world_map.demes["id"])[0][0]]
            d1_label = current_deme_ids[np.where(connection["deme_1"] == world_map.demes["id"])[0][0]]
            if (d0_label is not None) and (d1_label is not None):
                output.write("\tp" + d0_label + ".setMigrationRates(p" + d1_label + ", " + str(default_migration_rate) + ");\n")
                output.write("\tp" + d1_label + ".setMigrationRates(p" + d0_label + ", " + str(default_migration_rate) + ");\n")
        output.write("}\n\n")
        current_suitabilities = suitabilities.copy()
        current_start_time = epoch
    
