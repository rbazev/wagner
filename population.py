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
        """
        Number of individual genotypes in the population.

        >>> founder = population.Population.generate_founder(5,.2,1)
        >>> g0 = population.Population.found_clonal_population(founder,10)
        >>> g0.population_size
        10
        """
        return len(self.organisms)      

    def set_population_size(self, population_size):
        self.population_size = population_size
        
    @staticmethod
    def found_clonal_population(founder,population_size):
        """
        Make a population by copying a founder genotype into a list population_size times.

        >>> founder = population.Population.generate_founder(5,.2,1)
        >>> g0 = population.Population.found_clonal_population(founder,10)
        >>> g0
        <population.Population at 0x344fb30>
        >>> g0.organisms[2].gene_network
        array([[ 0.48238438,  0.        ,  0.        ,  0.08485553,  0.49141507],
        [-0.79636121,  0.        ,  0.        ,  0.        ,  0.        ],
        [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ],
        [ 0.        ,  0.        ,  0.        ,  0.        , -0.88783723],
        [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ]])
        """
        new_population = Population()
        while new_population.population_size < population_size:
            new_population.organisms.append(copy.deepcopy(founder))
        return new_population

    @staticmethod
    def generate_founder(n_genes, connectivity):
        """
        Generate a founder consisting of a gene network with random interaction strengths and a mutation rate.

        >>> founder = population.Population.generate_founder(5,.2)
        >>> founder.gene_network
        array([[ 0.48238438,  0.        ,  0.        ,  0.08485553,  0.49141507],
        [-0.79636121,  0.        ,  0.        ,  0.        ,  0.        ],
        [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ],
        [ 0.        ,  0.        ,  0.        ,  0.        , -0.88783723],
        [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ]])
        """
        founder = genotype.Genotype.generate_random(n_genes, connectivity)
        founder.generate_random_initial_expression_state()
        return founder
    
    def get_population_fitness(self):
        '''
        Calculate fitness of each genotype, create a vector of each genotype's fitness in the population and make it an attribute of the population
        '''
        self.all_fitness = []
        for i in range(len(self.organisms)):
            self.organisms[i].develop()
            self.organisms[i].calculate_fitness()
            self.all_fitness.append(self.organisms[i].fitness) 
            
#    @staticmethod
#    def recombine (genotype1, genotype2):
#        newmatrix = np.zeros((genotype1.n_genes, genotype1.n_genes)) 
#        for x in range (0, genotype1.n_genes):        
#            if random.random() >= .5:  
#                newmatrix[x] = genotype1.gene_network[x]
#            else:
#                newmatrix[x] = genotype2.gene_network[x]
#        return Genotype(newmatrix) 

                        
    def choose_genotype(self): 
        '''
        Generates all offspring for the next generation, randomly, weighted by fitness
        '''
        cum_fitness = np.cumsum(self.all_fitness)
        rand_array = rnd.random_sample(self.population_size)
        mult_rand_array = np.multiply(rand_array, cum_fitness[self.population_size-1])
        indices = np.searchsorted(cum_fitness, mult_rand_array)
        chosen_genotypes = []
        for i in range(self.population_size):
            chosen_genotypes.append(self.organisms[indices[i]])
        return chosen_genotypes
        

    def asexually_reproduce_pop (self):
        '''
        Asexually reproduces the population with a probability of mutations
        '''
        new_pop = Population()
        new_pop.organisms = Population.choose_genotype(self)
        for i in range(len(new_pop.organisms)):
            new_pop.organisms[i].mutate_random(rnd.poisson(new_pop.organisms[i].mutation_rate))
        return new_pop
            
       
       
    def sexually_reproduce_pop (self): 
        '''
        Sexually reproduces the population, choosing 2*population parents, recombining their gene networks, and mutating the product
        '''
        new_pop = Population()
        parents1 = Population.choose_genotype(self)
        parents2 = Population.choose_genotype(self)
        for i in range(len(parents1)):
            recombined = Genotype.recombine(parents1[i], parents2[i])
            recombined.mutate_random(rnd.poisson(recombined.mutation_rate))
            new_pop.organisms.append(recombined)
        return new_pop
            
        

    

