"""
genotype.py: implementation of the gene network model described in Siegal & Bergman (2002) Waddington's canalization
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


__author__ = 'Ricardo Azevedo, Christina Burch, Kayla Peck, Amanda Whitlock'
__version__ = "0.0.1"

class Population(object):
    """ 
    A population of gene networks
    """
  
  def set_population_size(self, population_size):
    self.population_size = population_size
  
  def set_stabilizing_selection(self, stabilizing_selection_strength):
    self.stabilizing_selection_strength = stabilizing_selection_strength
        
    
  def sexually_reproduce_pop (self, N):
    new_pop = []
    organism = 0
    while organism < N:
      genotype1 = choose_genotype #choose_genotype is the random-weighted selection method
      genotype2 = choose_genotype
      rec_offspring = recombine(genotype1, genotype2)
      mut_offspring = rec_offspring.mutate_random(rnd.poisson(rec_offspring.mutation_rate))
      new_pop.append(mut_offspring)
      organism = organism + 1
    return new_pop
        
    
  def choose_genotype(self, genotype.fitness)
    #randomly choose the genotype to replicate, weighted by fitness
