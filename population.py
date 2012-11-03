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
    Attributes: population_size
    """
  @staticmethod
  def found_population(population_size):
    founding_pop = Population()
    founding_pop.population_size = population_size
    founder = Genotype.generate_random(Genotype.n_genes, Genotype.connectivity)
    founding_pop.organisms = []
    popsize = 0
    while organism_num < population_size:
      founding_pop.organisms.append(founder)
      organism_num = organism_num + 1
    return founding.pop 
  
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
        
   def choose_offspring (self):   
    sumfitness = 0      #sum of the fitnesses of every genotype in pop, to be used as top range for generating rand number
    for genotype in self.organisms:
      sumfitness += genotype.fitness     
    rand = random.randrange(0, sumfitness)
    total = 0       #will iterate over genotypes, summing fitnesses until reaching random number
    for genotype in self.organisms:
        total += Genotype.fitness
        if total > rand
          return genotype
    
  def get_next_generation_asexual(self):
    next_generation = Population() 
    next_generation.population_size = self.population_size
    next_generation_stabilizing_selection_strength = self.stabilizing_selection_strength
    next_generation.organisms = []
    organism_num = 0
    while organism_num < next_generation.population_size
      randchoice = self.choose_organism
      mut_offspring = Genotype.randchoice.mutate_random(rnd.poisson(Genotype.randchoice.mutation_rate))
      next_generation.organisms.append(mut_offspring)
      organism_num = organism_num + 1
    
  def get_next_generation_sexual(self):
    next_generation = Population() 
    next_generation.population_size = self.population_size
    next_generation_stabilizing_selection_strength = self.stabilizing_selection_strength
    next_generation.organisms = []
    organism_num = 0
    while organism_num < next_generation.population_size
      parent1 = self.choose_organism
      parent2 = self.choose_organism
      offspring = recombine (parent1, parent2)
      mut_offspring = Genotype.rec_offspring.mutate_random(rnd.poisson(Genotype.rec_offspring.mutation_rate))
      next_generation.organisms.append(mut_offspring)
      organism_num = organism_num + 1
