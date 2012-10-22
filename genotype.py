"""
genotype.py: implementation of the gene network model described in Siegal & Bergman (2002) Waddington's canalization
    revisited: developmental stability and evolution. PNAS 99: 10528-32.

Contains the Genotype class.

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

class Genotype(object):
    """ 
    A gene network defined as described in: Wagner (1996) Does evolutionary plasticity evolve? *Evolution* 50: 1008-1023.

    Attributes:

        adj_matrix: adjacency matrix representation of the gene network
        connectivity: connectivity density
        gene_network: matrix representation of the gene network such that columns represent regulators and rows represent
            their targets
        graph: PyGraphviz representation of the gene network
        n_genes: number of genes
        n_interactions: number of (nonzero) interactions in the gene network
        mean_abs_strength: mean absolute strength of interactions (excluding zeros)

    """

    def __init__(self, matrix):
        """
        Create Genotype object based on a matrix.  Check that the matrix is square and a numpy.ndarray.

        Parameters:

            matrix: matrix representation of the gene network such that columns represent regulators and rows represent
                their targets (type: numpy.ndarray)

        >>> net = Genotype(np.array([[-1.,  0., -1.], [ 0.,  3.,  0.], [ 0.,  2., 1.]]))
        >>> net.gene_network
        array([[-1.,  0., -1.],
               [ 0.,  3.,  0.],
               [ 0.,  2.,  1.]])
        """
        assert type(matrix) is np.ndarray
        assert matrix.shape[0] == matrix.shape[1]
        self.gene_network = matrix

    @property
    def n_genes(self):
        """
        Number of genes in the network.

        >>> net = Genotype.generate_random(10, .35)
        >>> net.n_genes
        10
        """
        return len(self.gene_network.diagonal())

    @property
    def graph(self):
        """NetworkX representation of the gene network."""
        g = nx.DiGraph(data = self.gene_network.transpose())
        return g

    @property
    def interactions(self):
        """
        List of interactions in the form (regulator, target).

        >>> net = Genotype(np.array([[-1.,  0., -1.], [ 0.,  3.,  0.], [ 0.,  2., 1.]]))
        >>> net.interactions
        [(0, 0), (1, 1), (1, 2), (2, 0), (2, 2)]
        """
        return self.graph.edges()

    @property
    def n_interactions(self):
        """
        Number of nonzero interactions between genes.

        >>> net = Genotype.generate_random(10, .42)
        >>> net.n_interactions
        42
        """
        return len(self.interactions)

    @property
    def connectivity(self):
        """
        Connectivity density of the network.

        >>> net = Genotype.generate_random(10, .55)
        >>> net.connectivity
        0.55
        """
        return  self.n_interactions / float(self.n_genes * self.n_genes)

#    @property
#    def mean_abs_strength(self):
#        """Mean absolute strength of interactions (excluding zeros)."""
#        return np.abs(self.gene_network).sum() / self.n_interactions
#
#    @property
#    def graph(self):
#        """PyGraphviz representation of the gene network"""
#        g = pgv.AGraph(directed = True, strict = False)
#        g.node_attr['fontname'] = 'helvetica'
#        g.node_attr['fontsize'] = 16
#        g.node_attr['shape'] = 'circle'
#        g.node_attr['color'] = 'gray'
#        g.node_attr['style'] = 'filled'
#        # add genes
#        for gene in range(self.n_genes):
#            g.add_node(gene)
#        # add interactions
#        for regulator in range(self.n_genes):
#            for target in range(self.n_genes):
#                weight = self.gene_network[target, regulator]
#                if weight:
#                    g.add_edge(str(regulator), str(target))
#                    e = g.get_edge(str(regulator), str(target))
#                    e.attr['penwidth'] = np.abs(weight) / self.mean_abs_strength
#                    e.attr['color'] = 'red'
#                    if weight < 0:
#                        e.attr['color'] = 'blue'
#                        e.attr['arrowhead'] = 'tee'
#        g.layout(prog = 'dot')
#        return g
#
#    def draw_graph(self, filename):
#        """Draw gene network using graphviz.  Output PNG file."""
#        self.graph.draw(filename)

    @property
    def connected_components(self):
        """
        Connected components.  Return list containing lists of nodes in each connected component.

        >>> net = Genotype(np.array([[-1.,  0., -1.], [ 0.,  3.,  0.], [ 0.,  0., 1.]]))
        >>> net.connected_components
        [[0, 2], [1]]
        """
        return nx.connected_components(self.graph.to_undirected())

    @property
    def is_connected(self):
        """
        Whether the gene network is connected (type: bool).

        >>> net = Genotype(np.array([[-1.,  0., -1.], [ 0.,  3.,  0.], [ 0.,  0., 1.]]))
        >>> net.is_connected
        False
        """
        return len(self.connected_components) == 1

    def draw_graph(self):
        """
        Draw gene network using matplotlib.

        Note:
            Genes are labelled from 0.  Thick ends touch target genes.
            Self interactions are shown as gray nodes.
        """
        node_cols = ['white'] * self.n_genes
        diag = self.gene_network.diagonal()
        for i in range(self.n_genes):
            if diag[i] != 0:
                node_cols[i] = 'gray'
        nx.draw(self.graph, node_color = node_cols)

    @staticmethod
    def generate_random(n_genes, connectivity):
        """
        Generate random genotype with a given n_genes and connectivity.  Sample weights from standard normal distribution.
        Ensure that connectivity is between 0 and 1.

        Note:
            The realized connectivity will be different from the entered value if n_genes * n_genes * connectivity is not
            an integer (see code example below).

        Parameters:
            n_genes: number of genes (type: int)
            connectivity: connectivity density (type: float)

        >>> net = Genotype.generate_random(4, .3)
        >>> net.connectivity
        0.25
        """
        assert 0 <= connectivity <= 1
        n_sites = n_genes * n_genes
        n_nonzero_sites = n_sites * connectivity
        flat_matrix = np.zeros(n_sites)
        flat_matrix[0:(n_nonzero_sites)] = rnd.normal(size = n_nonzero_sites)
        rnd.shuffle(flat_matrix)
        return Genotype(np.reshape(flat_matrix, (n_genes, n_genes)))

    def mutate(self, regulator, target, new_strength):
        """
        Replaces interaction strength by new value.

        Parameters:

            regulator: number of regulator gene
            target: number of target gene
            new_strength: new interaction strength

        >>> net = Genotype(np.array([[-1.,  0., -1.], [ 0.,  3.,  0.], [ 0.,  0., 1.]]))
        >>> net.mutate(1, 2, 9.)
        >>> net.gene_network
        array([[-1.,  0., -1.],
               [ 0.,  3.,  0.],
               [ 0.,  9.,  1.]])
        """
        self.gene_network[target, regulator] = new_strength

    def mutate_random(self, n_mutations):
        """
        Mutate a given number of random interactions.  Sample new interaction strengths from standard normal distribution.
        Ensure that (i) mutations can only affect nonzero interactions, and (ii) the same interaction cannot be changed
        twice.

        Parameters:

            n_mutations: number of mutations

        >>> net = Genotype.generate_random(4, .3)
        >>> net.mutate_random(5)
        >>> net.connectivity
        0.25
        """
        if n_mutations > self.n_interactions:
            n_mutations = self.n_interactions
        new_strengths = rnd.normal(size = n_mutations)
        mut_interactions = []
        while len(mut_interactions) < n_mutations:
            x = rnd.randint(0, self.n_interactions)
            if not x in mut_interactions:
                mut_interactions.append(x)
        for i in range(n_mutations):
            interaction = self.interactions[mut_interactions[i]]
            self.mutate(interaction[0], interaction[1], new_strengths[i])

    def set_mutation_rate(self, mutation_rate):
        """
        Set mutation rate per genome per generation, that is, the expected number of mutations per generation.

        >>> net = Genotype.generate_random(4, .33)
        >>> net.set_mutation_rate(1.2)
        >>> net.mutation_rate
        1.2
        """
        if mutation_rate > self.n_interactions:
            self.mutation_rate = self.n_interactions
        else:
            self.mutation_rate = mutation_rate
            
    def generate_asexual_offspring(self):
        """
        Generate copy of a Genotype, allowing mutations to occur.
        Return a new Genotype object.

        >>> net = Genotype.generate_random(4, .5)
        >>> net.set_mutation_rate(1.2)
        >>> daughter_net = net.generate_asexual_offspring()
        >>> daughter_net.connectivity
        0.5
        """
        offspring = copy.deepcopy(self)
        offspring.mutate_random(rnd.poisson(offspring.mutation_rate))
        return offspring
   
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

        
if __name__ == "__main__":
    import doctest
    doctest.testmod()
