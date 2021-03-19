import networkx as nx 
import matplotlib.pyplot as plt
import math
import random

# Global constants
alpha = 0.5
beta = 1 - alpha
cycles = 10
initial_pheromone = 1.0
evaporation_rate = 0.7 # ρ

def decode_idx_to_int(edge):
    x, y = edge.split('_')
    
    x = int(x)
    y = int(y)

    return x, y

def calculate_distance(i, j):
    x1, y1 = decode_idx_to_int(i)
    x2, y2 = decode_idx_to_int(j)

    distance = math.sqrt( (x1 - x2) ** 2 + (y1 - y2) ** 2) 
    print("distance i j: " + str(distance) + " " + str(i) + " " + str(j))
    return distance


# Calculate separately the values to compute the proability
def tao_eta_ij(i, j): 
    try:
        pheromone_ij = MainGraph[i][j]['weight'] # τ_ij: pheromone between node i and j - tau
    except:
        MainGraph.add_edge(i, j, weight=initial_pheromone)
        SupportGraph.add_edge(i, j, weight=0.0)
        pheromone_ij = 1.0
    pheromone_ij = pow(pheromone_ij, alpha)

    x1, y1 = decode_idx_to_int(i)
    x2, y2 = decode_idx_to_int(j)

    # η_ij: 1 / distance of i to j - etal
    desirability_ij = 1 / calculate_distance(i, j)
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
MainGraph = nx.DiGraph() # Main graph
SupportGraph = nx.DiGraph() # Support graph
DirectedGraph = nx.DiGraph() # Graph with only directioned links

# Construct a graph based in positions of Job x Machine, like the Table 1 of the article
MainGraph.add_node("0_0", makespan=0, visited=False)
SupportGraph.add_node("0_0", makespan=0, visited=False)
DirectedGraph.add_node("0_0", makespan=0, visited=False)
ijb = len(info_job_machine)
# F_idx = str(ijb + 1) + "_0"
# MainGraph.add_node(F_idx, makespan=0, visited=False)
# SupportGraph.add_node(F_idx, makespan=0, visited=False)
# DirectedGraph.add_node(F_idx, makespan=0, visited=False)

# Store the makespan and visited info in each node
for i in range(len(info_job_machine)):
    for j in range(len(info_job_machine[0])):
        operation_tag = str(i + 1) + "_" + str(j + 1)
        time = info_job_machine[i][j][1]
        time = int(time)
        # print("time of ij postion is: " + str(time))
        MainGraph.add_node(operation_tag, makespan=time, visited=False)
        SupportGraph.add_node(operation_tag, makespan=0, visited=False)
        DirectedGraph.add_node(operation_tag, makespan=0, visited=False)

# Create the connections (edges) between i to j
for j in range(n):
    MainGraph.add_edge("0_0", str(j + 1) + "_1", weight=initial_pheromone)
    SupportGraph.add_edge("0_0", str(j + 1) + "_1", weight=0.0)
    DirectedGraph.add_edge("0_0", str(j + 1) + "_1", weight=0.0)

    for i in range(m - 1):
        MainGraph.add_edge(str(j + 1) + "_" + str(i + 1), str(j + 1) + "_" + str(i + 2), weight=initial_pheromone)
        SupportGraph.add_edge(str(j + 1) + "_" + str(i + 1), str(j + 1) + "_" + str(i + 2), weight=0.0)
        DirectedGraph.add_edge(str(j + 1) + "_" + str(i + 1), str(j + 1) + "_" + str(i + 2), weight=0.0)

    # MainGraph.add_edge(str(j + 1) + "_" + str(m), F_idx, weight=initial_pheromone)
    # SupportGraph.add_edge(str(j + 1) + "_" + str(m), F_idx, weight=0.0)
    # DirectedGraph.add_edge(str(j + 1) + "_" + str(m), F_idx, weight=0.0)

# for i in range(m):
#     for j in range(n - 1):
#         MainGraph.add_edge(str(j + 1) + "_" + str(i + 1), str(j + 2) + "_" + str(i + 1), weight=initial_pheromone)
#         SupportGraph.add_edge(str(j + 1) + "_" + str(i + 1), str(j + 2) + "_" + str(i + 1), weight=0.0)


# Create the connections (edges) between conflicted machines 
for jb in conflicted_job_machine:
    for i in range(len(jb) - 1):
        for j in range(i + 1, len(jb)):
            MainGraph.add_edge(jb[i], jb[j], weight=initial_pheromone)
            MainGraph.add_edge(jb[j], jb[i], weight=initial_pheromone)
            SupportGraph.add_edge(jb[i], jb[j], weight=0.0)
            SupportGraph.add_edge(jb[j], jb[i], weight=0.0)

# ---------- Starting the solutions ---------- #
ants = n // 2
best_solution = [0] * cycles
# auxGraph = nx.DiGraph(MainGraph)

for cycle in range(cycles):
    ants_solution = [0] * ants # L_k length of the makespan for each ant

    for a in range(ants):
        # AntGraph = auxGraph.copy()
        # AntGraph.clear_edges()

        # Setting the parameters for the first node
        actual_node = "0_0"
        MainGraph.nodes[actual_node]['visited'] = True

        # Verify if exist a non-visited neighbor
        neighbors_of_actual_node = list(DirectedGraph.adj[actual_node]) # Verify from directed ONLY graph
        candidate_nodes_to_visit = [] # tabu nodes or nodes members of S
        visited_nodes = []
        visited_nodes.append(actual_node)

        for n in neighbors_of_actual_node:
            if MainGraph.nodes[n]['visited'] == False:
                candidate_nodes_to_visit.append(n)
        
        len_candidate = len(candidate_nodes_to_visit)

        # While the ant have a path to walk, continue the algorithm
        while (len_candidate > 0):
            # Calculate the prob for each accesbile neighbor node by the function
            probability_for_each_node = [0.0] * len_candidate
            sum_tabu = 0.0

            print("actual node: " + str(actual_node))
            print("candidate nodes: ")
            print(candidate_nodes_to_visit)

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
            # SupportGraph.nodes[next_node]['makespan'] = MainGraph.nodes[next_node]['makespan'] + SupportGraph.nodes[actual_node]['makespan']
            # SupportGraph[actual_node][next_node]['weight'] += (1 / SupportGraph.nodes[next_node]['makespan']) # 1 / cumulative makespan
            # print("Previous G node: " + str(MainGraph.nodes[next_node]['makespan']))
            # print("Actual S graph cumulative makespan: " + str(SupportGraph.nodes[next_node]['makespan']))
            # print("Actual S graph cumulative pheremone: " + str(SupportGraph[actual_node][next_node]['weight']))

            # Store the makespan (solution) acquired for each ant
            # ants_solution[a] = max(ants_solution[a], SupportGraph.nodes[next_node]['makespan'])

            # Update the initial variables for next step to find the forward path
            actual_node = next_node
            MainGraph.nodes[actual_node]['visited'] = True

            neighbors_of_actual_node = list(DirectedGraph.adj[actual_node])
            # I need to remove used node in this loop
            candidate_nodes_to_visit.remove(actual_node) # tabu nodes
            visited_nodes.append(actual_node)

            for n in neighbors_of_actual_node:
                if MainGraph.nodes[n]['visited'] == False:
                    candidate_nodes_to_visit.append(n)
            
            # print(candidate_nodes_to_visit)
            len_candidate = len(candidate_nodes_to_visit)

        print("Order for visited nodes: " + str(visited_nodes))


        # for v in visited_nodes:
        #     # print(v)
        #     auxGraph.nodes[v]['visited'] = True

        #     if v == "0_0":
        #         AntGraph.add_edge("0_0", visited_nodes[1])
        #     else:
        #         # List non-visited neighbors
        #         neighbors_of_actual_node = list(auxGraph.adj[v]) # Verify from directed ONLY graph
        #         print(str(v) + ": " + str(neighbors_of_actual_node))
        #         for n in neighbors_of_actual_node:
        #             if auxGraph.nodes[n]['visited'] == False:
        #                 # SupportGraph.nodes[n]['visited'] = True
        #                 AntGraph.add_edge(v, n)

        # long_ant = nx.dag_longest_path(AntGraph)
        # print(long_ant)
        # a = 0
        # for l in long_ant:
        #     print(l)
        #     a += MainGraph.nodes[l]['makespan']

        # print(a)

        # After ONE ant visit the graph, reset the visited nodes in graph G
        for node_label, visited in MainGraph.nodes(data="visited"):
            if visited is not None:
                MainGraph.nodes[node_label]['visited'] = False

        # Creating a new graph in topological sort from visited nodes
        # for node_label, visited in AntGraph.nodes(data="visited"):
        #     if visited_nodes is not None:
        #         if node_label == "0_0":
        #             AntGraph.add_edge("0_0", visited_nodes[1])

    
        # After all ants visited the graph, reset the cumulative makespan of graph S
        for node_label, makespan in SupportGraph.nodes(data="makespan"):
            if makespan is not None:
                SupportGraph.nodes[node_label]['makespan'] = 0

    # Update pheromone of graph G
    for u, v, weight in MainGraph.edges(data="weight"):
        if weight is not None:
            MainGraph[u][v]['weight'] = MainGraph[u][v]['weight'] * evaporation_rate + SupportGraph[u][v]['weight']
            # print("Graph S pheromone: " + str(SupportGraph[u][v]['weight']))
            # print("Updated pheromone[" + str(u) + "][" + str(v) + "]: " + str(MainGraph[u][v]['weight']))
            # print("Graph S cumulative pheromone[" + str(u) + "][" + str(v) + "]: " + str(SupportGraph[u][v]['weight']))

    # After all ants visited the graph, reset the cumulative pheromone of graph S
    for u, v, weight in SupportGraph.edges(data="weight"):
        if weight is not None:
            # SupportGraph[u][v]['weight'] = 0.0
            SupportGraph[u][v]['weight'] = MainGraph[u][v]['weight']


# print(min(ants_solution))

# print(MainGraph.nodes())
# print(MainGraph.edges())
plt.subplot(221)
nx.draw_circular(MainGraph, with_labels=True)
plt.subplot(222)
nx.draw_circular(SupportGraph, with_labels=True)
plt.subplot(223)
nx.draw_circular(DirectedGraph, with_labels=True)
# plt.subplot(224)
# nx.draw_circular(AntGraph, with_labels=True)
plt.show()