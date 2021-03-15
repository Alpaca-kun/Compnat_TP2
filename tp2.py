import networkx as nx 
import matplotlib.pyplot as plt
import math
import random

# Global constants
alpha = 0.2
beta = 1 - alpha
cycles = 100
evaporation_rate = 0.7 # ρ

def beauty_print(arr):
    for i in arr:
        print(i)

def tao_eta_ij(i, j): # actual node and candidate node
    pheromone_ij = G[i][j]['weight'] # τ_ij: pheromone between node i and j - tau
    pheromone_ij = pow(pheromone_ij, alpha)

    x1, y1 = decode_idx_to_int(i)
    x2, y2 = decode_idx_to_int(j)

    # η_ij: 1 / distance of i to j - etal
    desirability_ij = 1 / math.sqrt( (x1 - x2) ** 2 + (y1 - y2) ** 2) 
    desirability_ij = pow(desirability_ij, beta)

    return pheromone_ij * desirability_ij

def amount_pheromone():
    pheromone_ij *= evaporation_rate # ρ * τ_ij(t - 1) 
    pheromone_k_ij = 1 / total_makespan # τ^k_ij 
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
        conflicted_job_machine[int(required_machine[i])].append(str(j + 1) + "_" + str(i + 1))

# ---------- Initializing the graph ----------
G = nx.DiGraph() # Main graph
S = nx.DiGraph() # Support graph

# Construct a graph based in positions of Job x Machine, like the Table 1 of the article
G.add_node("0_0", makespan=0, visited=False) # previous I node
S.add_node("0_0", makespan=0, visited=False)
ijb = len(info_job_machine)
F_idx = str(ijb + 1) + "_0"
G.add_node(F_idx, makespan=0, visited=False) # previous F node
S.add_node(F_idx, makespan=0, visited=False)

for i in range(len(info_job_machine)):
    for j in range(len(info_job_machine[0])):
        operation_tag = str(i + 1) + "_" + str(j + 1)
        time = info_job_machine[i][j][1]
        time = int(time)
        # print("time of ij postion is: " + str(time))
        G.add_node(operation_tag, makespan=time, visited=False)
        S.add_node(operation_tag, makespan=0, visited=False)
        # print(G.nodes[operation_tag])

for j in range(n):
    G.add_edge("0_0", str(j + 1) + "_1", weight=1.0)
    S.add_edge("0_0", str(j + 1) + "_1", weight=0.0)

    for i in range(m - 1):
        G.add_edge(str(j + 1) + "_" + str(i + 1), str(j + 1) + "_" + str(i + 2), weight=1.0)
        S.add_edge(str(j + 1) + "_" + str(i + 1), str(j + 1) + "_" + str(i + 2), weight=0.0)
    # G.add_edge(str(j) + "_" + str(m - 1), "F", weight=1.0)
    G.add_edge(str(j + 1) + "_" + str(m), F_idx, weight=1.0)
    S.add_edge(str(j + 1) + "_" + str(m), F_idx, weight=0.0)

for jb in conflicted_job_machine:
    for i in range(len(jb) - 1):
        for j in range(i + 1, len(jb)):
            G.add_edge(jb[i], jb[j], weight=1.0)
            G.add_edge(jb[j], jb[i], weight=1.0)
            S.add_edge(jb[i], jb[j], weight=0.0)
            S.add_edge(jb[j], jb[i], weight=0.0)

# ---------- Starting the solutions ---------- #
ants = n // 2
ants_solution = [] * ants # L_k length of the makespan for each ant

t = 0
while (t < cycles):
    for a in range(ants):
        actual_node = "0_0"
        G.nodes[actual_node]['visited'] = True
        number_of_neighbors = G.degree[actual_node]

        # while (number_of_neighbors > 0)
        neighbors_of_actual_node = list(G.adj[actual_node])
        # print("neighbors of I: " + str(neighbors_of_actual_node))
        
        # Verify the nodes already have been visited
        candidate_nodes = [] # tabu nodes
        for n in neighbors_of_actual_node:
            if G.nodes[n]['visited'] == False:
                candidate_nodes.append(n)

        # Calculate the prob for each node by the function
        len_candidate = len(candidate_nodes)
        probability_for_each_node = [0.0] * len_candidate
        sum_tabu = 0.0

        for c in range(len_candidate):                
            sum_tabu += tao_eta_ij(actual_node, candidate_nodes[c])

        for c in range(len_candidate):
            probability_for_each_node[c] = tao_eta_ij(actual_node, candidate_nodes[c]) / sum_tabu
            # print(probability_for_each_node[c])
    
        # sum_probability = sum(probability_for_each_node)
        # print("Sum is: " + str(sum_probability))

        # Roulette selection
        random_choice = random.random() # Random float:  0.0 <= x < 1.0
        cumulative_sum = [0.0] * len_candidate
        cumulative_sum[0] = probability_for_each_node[0]
        for c in range(len_candidate):
            if c == 0:
                pass
            cumulative_sum[c] = probability_for_each_node[c] + cumulative_sum[c - 1]
            # print("cumulative sum is :" + str(cumulative_sum[c]))

        selected_node_idx = 0
        for c in range(len_candidate):
            if c == 0:
                if random_choice <= cumulative_sum[c] and random_choice >= 0:
                    selected_node_idx = 0
                    break
            elif random_choice <= cumulative_sum[c] and random_choice > cumulative_sum[c - 1]:
                selected_node_idx = c
                break
        
        # print("random choice is: " + str(random_choice))
        # print("selected node is: " + str(selected_node_idx))

        # move the ant to next selected node
        next_node = candidate_nodes[selected_node_idx]

        # update makespan and pheromone in graph S
        S.nodes[next_node]['makespan'] = G.nodes[next_node]['makespan'] + G.nodes[actual_node]['makespan']
        S[actual_node][next_node]['weight'] += 1 / S.nodes[next_node]['makespan']
        # print(S.nodes[next_node]['makespan'])
        print(S[actual_node][next_node]['weight'])

        # Reset graph S
        for u, v, weight in S.edges(data="weight"):
            if weight is not None:
                # print(u, v, weight)
                # G[u][v]['weight'] *= evaporation_rate
                S[u][v]['weight'] = 0.0
                # print(A[u][v]['weight'])
                # print(G[u][v]['weight'])
                pass


    t += 1

# G['0_0']['1_1']['weight'] = 100
# print(G['0_0']['0_1']['weight']) # Start, end, weight between start and end

# print(G.nodes())
# print(G.edges())
plt.subplot(221)
nx.draw_circular(G, with_labels=True)
plt.subplot(222)
nx.draw_circular(S, with_labels=True)
plt.show()
