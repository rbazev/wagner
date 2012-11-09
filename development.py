"""
development.py

contains the Development class
"""    
import copy
import numpy as np
import numpy.random as rnd
import random
#import pygraphviz as pgv
import networkx as nx
import math
import genotype
import population

class Development(object):
    
    def calc_stab_sel_fitness(self, stabilizing_selection_strength):
        D = calculate_equilibrium_steady_state(average_expression_pattern, gene_expression_pattern, Genotype.n_genes) 
        self.fitness = math.log(D/stabilizing_selection_strength)    

    def generate_random_initial_expression_state(self):
        '''
        Generate an initial expression state - an array of size n_genes, filled randomly with 1 or -1
        '''
        self.initial_expression_state = np.round(rnd.random(self.n_genes))*2 - 1

    def set_activation_constant(self, activation_constant):
        '''
        Sets the activation constant and makes sure it is non-negative
        '''
        assert activation_constant > 0
        self.activation_constant = activation_constant

    @staticmethod
    def sigmoidal_filter_function(activation_constant, current_expression_state_index):
        '''
        Calculates the expression level of a gene by filtering the total regulatory input for the gene using the sigmoidal function
        f(x) = 2/(1+e^-ax) - 1. See Supplementary Figure 1.
        '''
        return (2/(1+math.exp(-activation_constant*current_expression_state_index)) - 1)

    def develop(self, n_steps):
        '''
        Simulates development - multiplies gene network R by initial expresssion state S(0) for n_steps number of times
        For each product, it is passed through the sigmoidal filter function (see supplementary information for 2006 paper)
        which then acts as S(t). The last 10 S states are assigned to the equilibrium expression state variable to be checked for stability.
        '''
        gene_expression_state = []
        gene_expression_state.append(self.initial_expression_state)
        for t in range(0,n_steps):
            current_expression_state = np.dot(self.gene_network, gene_expression_state[t])
            filtered_expression_state = []
            for x in range(0, self.n_genes):
                filtered_expression_state.append(Genotype.sigmoidal_filter_function(self.activation_constant,current_expression_state[x]))
            gene_expression_state.append(filtered_expression_state)
        self.gene_expression_pattern = np.array(gene_expression_state)

    def set_tau(self, tau):
        '''
        Sets the value of tau, the number of iterations that are included in calculating the equilibrium steady state
        '''
        assert tau > 0
        assert tau < len(self.gene_expression_pattern)
        self.tau = tau

    def average_expression_pattern(self):
        '''
        Calculates the average expression pattern from the last tau expression states
        '''
        interval_expression_pattern = self.gene_expression_pattern[(len(self.gene_expression_pattern)-self.tau):len(self.gene_expression_pattern)]
        self.average_expression_pattern = np.mean(interval_expression_pattern, axis=0)

    @staticmethod
    def calculate_equilibrium_steady_state(average_expression_pattern, gene_expression_pattern, n_genes):
        '''
        Calculates the value for the equilibrium steady state, to be used to check against the appropriate criterion (i.e. < 10^-3)
        '''
        difference_from_average_expression = (np.subtract(average_expression_pattern, gene_expression_pattern)**2)/(4*n_genes)
        return difference_from_average_expression

    @property
    def developmentally_stable(self):
        '''
        Checks the equilibrium expression state of the gene network. If the final sum is less than 10-3, it is stable.
        '''
        equilibrium_steady_state = []
        for x in range ((len(self.gene_expression_pattern) - self.tau),len(self.gene_expression_pattern)):
            equilibrium_steady_state.append(Genotype.calculate_equilibrium_steady_state(self.average_expression_pattern, self.gene_expression_pattern[x], self.n_genes))

        if np.sum(equilibrium_steady_state) < 0.001:
            return True
        else:
            return False

    def set_stabilizing_selection(self, stabilizing_selection_strength):
        self.stabilizing_selection_strength = stabilizing_selection_strength
            
    def calc_fitness(self, stabilizing_selection_strength):
        D = Development.calculate_equilibrium_steady_state(Development.average_expression_pattern, Development.gene_expression_pattern, n_genes)
        self.fitness = math.log(D/stabilizing_selection_strength)
        
        
        
        
        
        
        
        
        
        
