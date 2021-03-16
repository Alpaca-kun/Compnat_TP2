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

def decode_idx_to_int(edge):
    x, y = edge.split('_')
    
    x = int(x)
    y = int(y)

    return x, y

# Calculate separately the values to compute the proability
def tao_eta_ij(i, j): 
    pheromone_ij = G[i][j]['weight'] # τ_ij: pheromone between node i and j - tau
    pheromone_ij = pow(pheromone_ij, alpha)

    x1, y1 = decode_idx_to_int(i)
    x2, y2 = decode_idx_to_int(j)

    # η_ij: 1 / distance of i to j - etal
    desirability_ij = 1 / math.sqrt( (x1 - x2) ** 2 + (y1 - y2) ** 2) 
    desirability_ij = pow(desirability_ij, beta)

    return pheromone_ij * desirability_ij

# def amount_pheromone():
#     pheromone_ij *= evaporation_rate # ρ * τ_ij(t - 1) 
#     pheromone_k_ij = 1 / total_makespan # τ^k_ij 
#     total_pheromone = pheromone_ij + pheromone_k_ij
#     
#     return total_pheromone

# n = no of jobs, m = no of machines
n, m = [int(x) for x in input().split()]

info_job_machine = []
conflicted_job_machine = [[] for i in range(m)]

# Reading the input
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

# Store the makespan and visited info in each node
for i in range(len(info_job_machine)):
    for j in range(len(info_job_machine[0])):
        operation_tag = str(i + 1) + "_" + str(j + 1)
        time = info_job_machine[i][j][1]
        time = int(time)
        # print("time of ij postion is: " + str(time))
        G.add_node(operation_tag, makespan=time, visited=False)
        S.add_node(operation_tag, makespan=0, visited=False)

# Create the connections (edges) between i to j
for j in range(n):
    G.add_edge("0_0", str(j + 1) + "_1", weight=1.0)
    S.add_edge("0_0", str(j + 1) + "_1", weight=0.0)

    for i in range(m - 1):
        G.add_edge(str(j + 1) + "_" + str(i + 1), str(j + 1) + "_" + str(i + 2), weight=1.0)
        S.add_edge(str(j + 1) + "_" + str(i + 1), str(j + 1) + "_" + str(i + 2), weight=0.0)

    G.add_edge(str(j + 1) + "_" + str(m), F_idx, weight=1.0)
    S.add_edge(str(j + 1) + "_" + str(m), F_idx, weight=0.0)

# Create the connections (edges) between conflicted machines 
for jb in conflicted_job_machine:
    for i in range(len(jb) - 1):
        for j in range(i + 1, len(jb)):
            G.add_edge(jb[i], jb[j], weight=1.0)
            G.add_edge(jb[j], jb[i], weight=1.0)
            S.add_edge(jb[i], jb[j], weight=0.0)
            S.add_edge(jb[j], jb[i], weight=0.0)

# ---------- Starting the solutions ---------- #
ants = n // 2
best_solution = [0] * cycles

for cycle in range(cycles):
    ants_solution = [0] * ants # L_k length of the makespan for each ant

    for a in range(ants):
        # Setting the parameters for the first node
        actual_node = "0_0"
        G.nodes[actual_node]['visited'] = True

        # Verify if exist a non-visited neighbor
        neighbors_of_actual_node = list(G.adj[actual_node])
        candidate_nodes_to_visit = [] # tabu nodes

        for n in neighbors_of_actual_node:
            if G.nodes[n]['visited'] == False:
                candidate_nodes_to_visit.append(n)
        
        len_candidate = len(candidate_nodes_to_visit)

        # While the ant have a path to walk, continue the algorithm
        while (len_candidate > 0):
            # Calculate the prob for each accesbile neighbor node by the function
            probability_for_each_node = [0.0] * len_candidate
            sum_tabu = 0.0

            for c in range(len_candidate):                
                sum_tabu += tao_eta_ij(actual_node, candidate_nodes_to_visit[c])

            for c in range(len_candidate):
                probability_for_each_node[c] = tao_eta_ij(actual_node, candidate_nodes_to_visit[c]) / sum_tabu
        
            # Cumulative Roulette Selection
            random_choice = random.random() # Random float:  0.0 <= x < 1.0
            cumulative_sum = [0.0] * len_candidate
            cumulative_sum[0] = probability_for_each_node[0]

            for c in range(len_candidate):
                if c == 0:
                    pass
                cumulative_sum[c] = probability_for_each_node[c] + cumulative_sum[c - 1]
                # print("cumulative sum is :" + str(cumulative_sum[c]))

            selected_node_idx_to_move = 0
            for c in range(len_candidate):
                if c == 0:
                    if random_choice <= cumulative_sum[c] and random_choice >= 0:
                        selected_node_idx_to_move = 0
                        break
                elif random_choice <= cumulative_sum[c] and random_choice > cumulative_sum[c - 1]:
                    selected_node_idx_to_move = c
                    break
            
            # print("random choice is: " + str(random_choice))
            # print("selected node is: " + str(selected_node_idx_to_move))

            # Move the ant to the next selected node
            next_node = candidate_nodes_to_visit[selected_node_idx_to_move]
            # print("actual node idx: " + str(actual_node))
            # print("next node idx: " + str(next_node))

            # Update makespan and pheromone in graph S
            S.nodes[next_node]['makespan'] = G.nodes[next_node]['makespan'] + S.nodes[actual_node]['makespan']
            S[actual_node][next_node]['weight'] += (1 / S.nodes[next_node]['makespan']) # 1 / cumulative makespan
            # print("Previous G node: " + str(G.nodes[next_node]['makespan']))
            # print("Actual S graph cumulative makespan: " + str(S.nodes[next_node]['makespan']))
            # print("Actual S graph cumulative pheremone: " + str(S[actual_node][next_node]['weight']))

            # Store the makespan (solution) acquired for each ant
            ants_solution[a] = max(ants_solution[a], S.nodes[next_node]['makespan'])

            # Update the initial variables for next step to find the forward path
            actual_node = next_node
            G.nodes[actual_node]['visited'] = True

            neighbors_of_actual_node = list(G.adj[actual_node])
            candidate_nodes_to_visit = [] # tabu nodes

            for n in neighbors_of_actual_node:
                if G.nodes[n]['visited'] == False:
                    candidate_nodes_to_visit.append(n)
            
            len_candidate = len(candidate_nodes_to_visit)

        # print(ants_solution)

        # After ONE ant visit the graph, reset the visited nodes in graph G
        for node_label, visited in G.nodes(data="visited"):
            if visited is not None:
                G.nodes[node_label]['visited'] = False
    
        # After all ants visited the graph, reset the cumulative makespan of graph S
        for node_label, makespan in S.nodes(data="makespan"):
            if makespan is not None:
                S.nodes[node_label]['makespan'] = 0

    # Update pheromone of graph G
    for u, v, weight in G.edges(data="weight"):
        if weight is not None:
            G[u][v]['weight'] = G[u][v]['weight'] * evaporation_rate + S[u][v]['weight']
            # print("Graph S pheromone: " + str(S[u][v]['weight']))
            # print("Updated pheromone[" + str(u) + "][" + str(v) + "]: " + str(G[u][v]['weight']))
            # print("Graph S cumulative pheromone[" + str(u) + "][" + str(v) + "]: " + str(S[u][v]['weight']))

    # After all ants visited the graph, reset the cumulative pheromone of graph S
    for u, v, weight in S.edges(data="weight"):
        if weight is not None:
            S[u][v]['weight'] = 0.0


print(min(ants_solution))

# G['0_0']['1_1']['weight'] = 100
# print(G['0_0']['0_1']['weight']) # Start, end, weight between start and end

# print(G.nodes())
# print(G.edges())
# plt.subplot(221)
# nx.draw_circular(G, with_labels=True)
# plt.subplot(222)
# nx.draw_circular(S, with_labels=True)
# plt.show()