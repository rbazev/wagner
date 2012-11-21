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
        genotype.Genotype.generate_random_initial_expression_state(founder)
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


    def set_percent_sexual (self, percent):
        '''
        Create a population with both sexual and asexual members
        '''
        assert 0 <= percent <=1
        self.percent_sexual = percent
        number_sexual = int(self.population_size*percent)
        for i in range(0, self.population_size):
            if i < number_sexual:
                self.organisms[i].sex_locus = "sexual"
            else:
                self.organisms[i].sex_locus = "asexual"
        
    def get_percent_sexual (self):
        ''' 
        Calculate the percentage of organisms that are sexual
        '''
        num_sexual = 0
        for i in range(0, self.population_size):
            if self.organisms[i].sex_locus == "sexual":
                num_sexual = num_sexual + 1
            else:
                next
        self.percent_sexual = num_sexual/self.population_size      
        
                            
    def choose_genotype(self): 
        '''
        Generates all offspring for the next generation, randomly, weighted by fitness
        '''
        cum_fitness = np.cumsum(self.all_fitness)
        rand_array = rnd.random_sample(self.population_size)
        mult_rand_array = np.multiply(rand_array, cum_fitness[self.population_size-1])
        indices = np.searchsorted(cum_fitness, mult_rand_array)
        chosen_genotypes = []
        for i in range(0,self.population_size):
            chosen_genotypes.append(self.organisms[indices[i]])
        return chosen_genotypes
        

    def reproduce_pop (self):
        '''
        Determines which reproduction method is appropriate.
        '''
        if self.percent_sexual == 1:
            new_pop = Population.sexually_reproduce_pop(self)
        elif self.percent_sexual == 0:
            new_pop = Population.asexually_reproduce_pop(self)
        else:
            new_pop = Population.reproduce_mixed_pop(self)
        return new_pop
            

    def reproduce_mixed_pop(self): 
        '''
        If the population contains both sexual and asexual organisms, new organisms are chosen and sex locus values compared.
        Sexual reproduction is dominant.
        '''
        new_pop = Population()
        chosen_genotypes1 = Population.choose_genotype(self)
        chosen_genotypes2 = Population.choose_genotype(self)
        for i in range(0,len(chosen_genotypes1)):
            if len(new_pop.organisms) < self.population_size:
                if (chosen_genotypes1[i].sex_locus == "asexual" and chosen_genotypes2[i].sex_locus == "asexual"):
                    chosen_genotypes1[i].mutate_random(rnd.poisson(chosen_genotypes1[i].mutation_rate))
                    new_pop.organisms.append(chosen_genotypes1[i])
                elif chosen_genotypes1[i].sex_locus == "sexual" and chosen_genotypes2[i].sex_locus == "sexual":
                    recombined = (genotype.Genotype.recombine(chosen_genotypes1[i], chosen_genotypes2[i]))
                    recombined.mutate_random(rnd.poisson(recombined.mutation_rate))
                    new_pop.organisms.append(recombined)
                    
                else:
                    recombined = (genotype.Genotype.recombine(chosen_genotypes1[i], chosen_genotypes2[i]))
                    recombined.mutate_random(rnd.poisson(recombined.mutation_rate))
                    new_pop.organisms.append(recombined)       
        return new_pop                    
     

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
            
        

    

