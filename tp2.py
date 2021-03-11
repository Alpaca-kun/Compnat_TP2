import networkx as nx 

def beauty_print(arr):
    for i in arr:
        print(i)

# n = no of jobs, m = no of machines
n, m = [int(x) for x in input().split()]

jobs = []
job_machine = [[] for i in range(m)]
for j in range(n):
    sequences = input().split()
    required_machine = sequences[0:][::2]
    processing_time = sequences[1:][::2]
    jobs.append([])

    for i in range(m):
        # job j has a list of required machine and processing time
        jobs[j].append((required_machine[i], processing_time[i]))
        job_machine[int(required_machine[i])].append(str(j) + "_" + str(i))

# beauty_print(job_machine)
print()
beauty_print(jobs)

# Initializing the graph
G = nx.DiGraph()
G.add_node("I")
G.add_node("F")

job_number = 0
for job in jobs:
    # print(job)
    seq_number = 0
    for s in job:
        operation_tag = str(job_number) + "_" + str(seq_number)
        G.add_node(operation_tag, time=s[1])
        # print(G.nodes[operation_tag])
        seq_number += 1
    job_number += 1
# print(G.nodes())
# print(G.edges())

for j in range(n):
    # print("I", str(j) + "_0")
    G.add_edge("I", str(j) + "_0")
    for i in range(m - 1):
        # print(str(j) + "_" + str(i), str(j) + "_" + str(i+1))
        G.add_edge(str(j) + "_" + str(i), str(j) + "_" + str(i+1))
    G.add_edge(str(j) + "_" + str(m - 1), "F")

# Fazer uma funcao de print bonito
# Montar as arestas dificeis
for jb in job_machine:
    for i in range(len(jb) - 1):
        for j in range(i + 1, len(jb)):
            # print(jb[i],jb[j])
            G.add_edge(jb[i], jb[j])
            G.add_edge(jb[j], jb[i])

# print(G.edges())