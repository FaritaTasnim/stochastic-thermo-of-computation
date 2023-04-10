import networkx as nx
import numpy as np

'''
Generate networks according to your modular hierarchy generative method.
Set N, a number of nodes, e.g., 100

Assign nodes to layers 0 - L-1 according to (a^L - 1)/(a - 1) = N, with a^l nodes in the lth layer.
Then within each layer, split apart the nodes into C_l communities. 

Then, assign edges according to p_up, p_down, p_intra, p_inter. 
Then use networkx to draw the network (potentially need to use a different software )
'''

a = 3
L = 5
N = (a**L - 1)/(a - 1)
nodes = np.arange(N)
layers = [nodes[int((a**(l) - 1)/(a - 1)):int((a**(l+1) - 1)/(a - 1))] for l in range(L)]

print (layers)
print([len(k) for k in layers])