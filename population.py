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
import development


__author__ = 'Ricardo Azevedo, Christina Burch, Kayla Peck, Amanda Whitlock'
__version__ = "0.0.1"

class Population(object):
    """
    A population of gene networks
    Attributes: population_size
    """
    def __init__(self):
        self.organisms = []

    @property
    def population_size(self):
        return len(self.organisms)
        
    def set_population_size(self, population_size):
        self.population_size = population_size

    @staticmethod
    def found_population(population_size):
        founding_pop = Population()
        founding_pop.population_size = population_size
        founding_pop.organisms = []
        while len(founding_pop.organisms) < population_size:
            founding_pop.organisms.append(copy.deepcopy(founder))
        return founding_pop

    @staticmethod
    def generate_founder(n_genes, connectivity, mutation_rate):
#        Genotype.n_genes = n_genes
#        Genotype.connectivity = connectivity
         founder = Genotype.generate_random(n_genes, connectivity)
#        founder.generate_random_initial_expression_state()
         founder.set_mutation_rate(mutation_rate)
#        founder.set_activation_constant(activation_constant)
#        founder.develop(dev_steps)
#        founder.set_tau(tau)
#        founder.developmentally_stable()
         return founder
    
    def get_fitness (self):
        self.all_fitness = []
        for i in range(len(founding_pop.organisms)):
            self.all_fitness.append(self.organisms[i].fitness) 
            
    @staticmethod
    def recombine(genotype1,genotype2):
        gene_num = genotype1.n_genes
        newmatrix = np.zeros((gene_num, gene_num)) 
        for x in range(0, gene_num): #iterate through loop * n_genes, ie for each row in the matrix
            parent = random.random()  #randomly pick which parent will contribute each row of the matrix
            if parent >= .5:
                chosen_genotype = genotype1
            else:
                chosen_genotype = genotype2
            row = chosen_genotype.gene_network[x] #set row x of the parent chosen equal to row
            newmatrix[x] = row #add row x of the parent chosen to the offspring's genotype
        offspring = Genotype(newmatrix)
        return offspring    
            
                
    def choose_genotype(self):
        cum_fitness = np.cumsum(self.all_fitness)
        rand_array = rnd.random_sample(self.population_size)
        mult_rand_array = np.multiply(rand_array, cum_fitness[self.population_size-1])
        indices = np.searchsorted(cum_fitness, mult_rand_array)
        chosen_genotypes = []
        for i in range(self.population_size):
            chosen_genotypes.append(self.organisms[indices[i]])
        return chosen_genotypes
        
        
    def sexually_reproduce_pop (self):
        new_pop = Population()
        new_pop.offspring = []
        while len(new_pop.offspring) < self.population_size:
            genotype1 = choose_genotype #choose_genotype is the random-weighted selection method
            genotype2 = choose_genotype
            rec_offspring = recombine(genotype1, genotype2)
            mut_offspring = rec_offspring.mutate_random(rnd.poisson(rec_offspring.mutation_rate))
            new_pop.append(mut_offspring)
        return new_pop
                
    
    def asexually_reproduce_pop (self): 
        self.organsisms = [] #this will make it faster but is there ever a time when we will want to preserve the population?
        new_pop = copy.deepcopy(self)
        chosen_genotypes = choose_genotype(self)
        for i in range(len(chosen_genotypes)):
            mut_offspring = chosen_genotypes[i].mutate_random(rnd.poisson(Genotype.chosen_genotypes[i].mutation_rate)) 
            new_pop.organisms(append.mut_offspring)
        return new_pop
        
       
    def sexually_reproduce_pop (self): 
        self.organsisms = [] #this will make it faster but is there ever a time when we will want to preserve the population?
        new_pop = copy.deepcopy(self)
        chosen_genotypes1 = choose_genotype(self)
        chosen_genotypes2 = choose_genotype(self)
        recombined_offspring = []
        for i in range(len(chosen_genotypes1)):
            rec_offsping = recombine(chosen_genotype1[i], chosen_genotype2[i])
            recombined_offspring.append.rec_offspring
        for i in range(len(chosen_genotypes2)):
            mut_offspring = recombined_offspring[i].mutate_random(rnd.poisson(Genotype.mut_offspring[i].mutation_rate))
            new_pop.organisms.append.mut_offspring
        return new_pop
        

    

