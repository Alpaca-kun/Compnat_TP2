import networkx as nx 
import matplotlib.pyplot as plt
import math
import random

# Global constants
alpha = 0.2
beta = 1 - alpha
cycles = 1000
evaporation_rate = 0.7 # ρ

def beauty_print(arr):
    for i in arr:
        print(i)

def probability(x1, x2, y1, y2):
    pheromone_ij# τ_ij: pheromone between node i and j - tau
    # η_ij: 1 / distance of i to j - etal
    desirability_ij = 1 / math.sqrt( (x1 - x2) ** 2 + (y1 - y2) ** 2) 

def amount_pheromone():
    pheromone_ij *= evaporation_rate # ρ * τ_ij(t - 1) 
    pheromone_k_ij = 1 / length_of_path # τ^k_ij 
    total_pheromone = pheromone_ij + pheromone_k_ij
    
    return total_pheromone

def decode_idx_to_int(edge):
    x, y = edge.split('_')
    
    x = int(x)
    y = int(y)

    return x, y

# n = no of jobs, m = no of machines
n, m = [int(x) for x in input().split()]

info_job_machine = []
conflicted_job_machine = [[] for i in range(m)]

for j in range(n):
    sequences = input().split()
    required_machine = sequences[0:][::2]
    processing_time = sequences[1:][::2]
    info_job_machine.append([])

    for i in range(m):
        # job j has a list of required machine and processing time
        info_job_machine[j].append((required_machine[i], processing_time[i]))
        # store the j-th job and i-th machine to verify the conflicted jobs that use the same machine
        conflicted_job_machine[int(required_machine[i])].append(str(j) + "_" + str(i))
        # print(conflicted_job_machine[int(required_machine[i])])

# ---------- Initializing the graph ----------
G = nx.DiGraph()
G.add_node("I", makespan=0)
G.add_node("F", makespan=0)
print("makespan is: " + int(G.nodes['I']['makespan']))

# Construct a graph based in positions of Job x Machine, like the Table 1 of the article
for i in range(len(info_job_machine)):
    for j in range(len(info_job_machine[0])):
        operation_tag = str(i) + "_" + str(j)
        time = info_job_machine[i][j][1]
        # print("time of ij postion is: " + str(time))
        G.add_node(operation_tag, makespan=time)
        # print(G.nodes[operation_tag])

for j in range(n):
    G.add_edge("I", str(j) + "_0", weight=1.0)
    for i in range(m - 1):
        G.add_edge(str(j) + "_" + str(i), str(j) + "_" + str(i+1), weight=1.0)
    G.add_edge(str(j) + "_" + str(m - 1), "F", weight=1.0)

for jb in conflicted_job_machine:
    for i in range(len(jb) - 1):
        for j in range(i + 1, len(jb)):
            G.add_edge(jb[i], jb[j], weight=1.0)
            G.add_edge(jb[j], jb[i], weight=1.0)

# ---------- Starting the solutions ---------- #
ants = n // 2
ants_solution = [] * ants # L_k length of the makespan for each ant

t = 0
while (t < cycles):
    for a in range(ants):
        # A = nx.DiGraph()

        actual_node = 'I'
        next_node = '' # empty string
        edge_selection = ''

        verify_neighbors = G.degree[actual_node]

        if t == 0:
            edges_of_I = G.degree['I']
            adj_of_I = list(G.adj['I'])
            # print("edges of I: " +str(edges_of_I))
            # print("neighbors of I: " +str(list(G.adj['I'])))
            
            # the prob for the Initial node is equivalent
            edge_selection = random.randrange(edges_of_I) # select the node
            next_node = str(adj_of_I[edge_selection])# store the next node 
            # print(next_node)
            # indexs = next_node.split('_')
            next_x, next_y = decode_idx_to_int(next_node)

            # A.add_edge('I', next_node)
            # nx.draw_circular(A, with_labels=True)

        else: 
            edge_selection = random.randrange(10)


    t += 1

# for g in G:
#   print(G[str(g)])

# for n, nbrsdict in G.adjacency():
#     for nbr, eattr in nbrsdict.items():
#         if "weight" in eattr:
#             # Do something useful with the edges
#             eattr['weight'] = 1000
#             print(eattr)

for u, v, weight in G.edges(data="weight"):
    if weight is not None:
        # print(u, v, weight)
        G[u][v]['weight'] *= evaporation_rate
        # print(G[u][v]['weight'])


G['0_0']['0_1']['weight'] = 100
print(G['0_0']['0_1']['weight']) # Start, end, weight between start and end

# print(G.nodes())
# print(G.edges())
# nx.draw_circular(G, with_labels=True)
# nx.draw_circular(A, with_labels=True)
# plt.show()
