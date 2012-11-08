"""
population.py: implementation of the gene network model described in Siegal & Bergman (2002) Waddington's canalization
    revisited: developmental stability and evolution. PNAS 99: 10528-32.

Contains the Population class.

"""

import copy
import numpy as np
import numpy.random as rnd
import random
#import pygraphviz as pgv
import networkx as nx
import math
import genotype

__author__ = 'Ricardo Azevedo, Christina Burch, Kayla Peck, Amanda Whitlock'
__version__ = "0.0.1"

class Population(object):
    """
    A population of gene networks
    Attributes: population_size
    """
    def __init__(self):
        self.population = []

    @property
    def population_size(self):
        return len(self.population)

    @staticmethod
    def found_population(population_size, n_genes, connectivity):
        founding_pop = Population()
        founding_pop.organisms = []
        popsize = 0
        organism_num = 0
        while organism_num < population_size:
            founding_pop.organisms.append(copy.deepcopy(Genotype.founder))
            organism_num = organism_num + 1
        return founding_pop

    @staticmethod
    def generate_founder(n_genes, connectivity, mutation_rate):
#        Genotype.n_genes = n_genes
#        Genotype.connectivity = connectivity
        founder = genotype.Genotype.generate_random(n_genes, connectivity)
#        founder.generate_random_initial_expression_state()
        founder.set_mutation_rate(mutation_rate)
#        founder.set_activation_constant(activation_constant)
#        founder.develop(dev_steps)
#        founder.set_tau(tau)
#        founder.developmentally_stable()
        return founder
  
    def set_population_size(self, population_size):
        self.population_size = population_size
  
    def set_stabilizing_selection(self, stabilizing_selection_strength):
        self.stabilizing_selection_strength = stabilizing_selection_strength
        
    
    def sexually_reproduce_pop (self, N):
        new_pop = Population()
        organism = 0
        while organism < N:
            genotype1 = choose_genotype #choose_genotype is the random-weighted selection method
            genotype2 = choose_genotype
            rec_offspring = recombine(genotype1, genotype2)
            mut_offspring = rec_offspring.mutate_random(rnd.poisson(rec_offspring.mutation_rate))
            new_pop.append(mut_offspring)
            organism = organism + 1
        return new_pop
        
    

