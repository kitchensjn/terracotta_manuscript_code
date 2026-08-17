import sys
sys.path.append("/Users/jameskitchens/Documents/GitHub/terracotta")
sys.path.append("../..")
import sim_utils
import pandas as pd
import matplotlib.pyplot as plt


sim_utils.create_trees_files(
    demes_path="dataset/demes.tsv",
    connections_path="dataset/connections.tsv",
    samples_path="dataset/samples.tsv",
    number_of_trees=10,
    pop_size=15,
    migration_rate=0.001,
    end_time=121000,
    output_directory="dataset"
)