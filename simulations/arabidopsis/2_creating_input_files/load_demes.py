#import demes
#import demesdraw

#tct = demes.load("demes_test.yaml")
#ax = demesdraw.tubes(tct)
#ax.figure.savefig("test.svg")


import fwdpy11
import numpy as np
from dataclasses import dataclass


@dataclass
class Parents:
    generation: int
    parents: np.array


@dataclass
class RecordParents:
    parents: list

    def __call__(self, pop, _):
        md = np.array(pop.diploid_metadata, copy=False)
        self.parents.append(Parents(pop.generation, np.copy(md["parents"])))

with open("demes_test.yaml", "r") as yaml:
    demography = fwdpy11.ForwardDemesGraph.from_demes(yaml.read(), burnin=10000, burnin_is_exact=True)
    pop = fwdpy11.DiploidPopulation(demography.initial_sizes, 1000.)
    rng = fwdpy11.GSLrng(54321)
    pdict = {'simlen': demography.final_generation,
         'demography': demography,
         # Everything below is standard
         # and covered in other vignettes
         'recregions': [],
         'nregions': [],
         'sregions': [],
         'rates': (0,0,0),
         'gvalue': fwdpy11.Multiplicative(2.)
         }
    params = fwdpy11.ModelParams(**pdict)
    recorder = RecordParents([])
    fwdpy11.evolvets(rng, pop, params, simplification_interval=1000000, recorder=recorder)
    
    #ts = pop.dump_tables_to_tskit()
    #print(ts.draw_text())
