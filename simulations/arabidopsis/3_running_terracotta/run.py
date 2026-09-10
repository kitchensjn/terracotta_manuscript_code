import sys
sys.path.append("/Users/jameskitchens/Documents/GitHub/terracotta")
import terracotta as tct


result, loglikelihood = tct.run(
    demes_path="../2_creating_input_files/dataset/demes.tsv",
    connections_path="../2_creating_input_files/dataset/connections.tsv",
    samples_path="../2_creating_input_files/dataset/samples.tsv",
    trees_dir_path="../2_creating_input_files/dataset/trees_small",
    output_file="output_no_coal.tsv"
)

print(result)
print(loglikelihood)