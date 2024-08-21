import itertools
import igraph as ig
import matplotlib.pyplot as plt
import copy
from collections import defaultdict, deque
class Graph:
    def __init__(self):
        self.nodes = []  # Each node is a list [node_id, value]
        self.edges = []  # Each edge is a list [[node1_id, node2_id], value]
    
    def add_node(self, node_id, value="C"):
        if node_id not in [node[0] for node in self.nodes]:
            self.nodes.append([node_id, value])
    
    def add_edge(self, node1_id, node2_id, value=1):
        self.edges.append([[node1_id, node2_id], value])

def terminal_carbon(graph):
    output = []
    for node in graph.nodes:
        node_id = node[0]
        # Count the number of edges connected to this node
        adjacent_count = sum(1 for edge in graph.edges if node_id in edge[0])
        
        # A terminal node has exactly one connection
        if adjacent_count == 1:
            output.append(node_id)
    
    return output

def all_chain_all(graph, start_node):
    if len(graph.nodes) == 1:
        return [[[graph.nodes[0][0]], []]]
    all_chain_output = []
    sub_chain_output = []
    def all_chain(chain, compound, sub_chain):
        end_reached = True
        atom_iter = []
        for atom in compound.nodes:
            atom_id = atom[0]
            if atom_id not in chain and {atom_id, chain[-1]} in [set(item[0]) for item in compound.edges]:
                atom_iter.append(atom_id)
        if len(atom_iter) == 0:
            all_chain_output.append(chain)
            sub_chain_output.append(sub_chain)
        elif len(atom_iter) == 1:
            atom_id = atom_iter[0]
            all_chain(chain + [atom_id], compound, sub_chain)
        else:
            for atom_id in atom_iter:
                all_chain(chain + [atom_id], compound, sub_chain + [len(chain)]*(len(atom_iter)-1))
    if start_node is None:
        for atom in terminal_carbon(graph):
            all_chain([atom], graph, [])
    else:
        all_chain([start_node], graph, [])
    #print("JI")
    return [list(x) for x in zip(all_chain_output, sub_chain_output)]
count_hydrogen = 6
def add_hydrogens(graph):
    global count_hydrogen
    hydrogens = []  # To keep track of hydrogen nodes added
    for node in graph.nodes:
        node_id, value = node
        if value == "C":
            # Count existing bonds
            existing_bonds = sum([edge[1] for edge in graph.edges if node_id in edge[0]])
            needed_hydrogens = 4 - existing_bonds
            
            for _ in range(needed_hydrogens):
                hydrogen_id = count_hydrogen
                count_hydrogen += 1
                graph.add_node(hydrogen_id, value="H")
                graph.add_edge(node_id, hydrogen_id)
                hydrogens.append(hydrogen_id)
    return hydrogens

def draw_graph(graph):
    g = ig.Graph(directed=False)
    
    # Add vertices with labels
    vertex_map = {}
    for node in graph.nodes:
        node_id, value = node
        v = g.add_vertex(name=node_id)
        vertex_map[node_id] = v
    
    # Add edges with corresponding colors
    edge_colors = []
    for edge in graph.edges:
        node1, node2 = edge[0]
        g.add_edge(vertex_map[node1], vertex_map[node2])
        edge_value = edge[1]
        
        if edge_value == 1:
            edge_colors.append("gray")
        elif edge_value == 2:
            edge_colors.append("red")
        elif edge_value == 3:
            edge_colors.append("green")
        else:
            edge_colors.append("black")  # Default color for unknown values
    
    # Set vertex labels
    vertex_labels = [node[1] for node in graph.nodes]
    g.vs['label'] = vertex_labels
    
    # Draw the graph
    layout = g.layout("fr")  # Fruchterman-Reingold layout
    ig.plot(g, layout=layout, vertex_label=g.vs['label'], vertex_size=30, vertex_color="lightblue", edge_color=edge_colors, edge_width=2, bbox=(800, 600), margin=50, target="graph_igraph.png")
    
    # Display the graph using matplotlib
    img = plt.imread('graph_igraph.png')
    plt.figure(figsize=(12, 8))
    plt.imshow(img)
    plt.axis('off')
    plt.show()

def explorable_subgraph(graph, start_node, forbidden_nodes):
    # Create a new Graph instance for the explorable subgraph
    subgraph = Graph()
    
    # Perform BFS to find the explorable subgraph
    if start_node in [node[0] for node in graph.nodes]:
        visited = set()
        queue = [start_node]
        
        while queue:
            node = queue.pop(0)
            if node not in visited and node not in forbidden_nodes:
                visited.add(node)
                subgraph.add_node(node, value=[n[1] for n in graph.nodes if n[0] == node][0])
                
                # Add edges to the subgraph if they connect to valid nodes
                for edge in graph.edges:
                    (n1, n2), value = edge
                    if n1 == node and n2 not in visited and n2 not in forbidden_nodes:
                        queue.append(n2)
                        if n2 not in forbidden_nodes:
                            subgraph.add_node(n2, value=[n[1] for n in graph.nodes if n[0] == n2][0])
                            subgraph.add_edge(n1, n2, value)
                    elif n2 == node and n1 not in visited and n1 not in forbidden_nodes:
                        queue.append(n1)
                        if n1 not in forbidden_nodes:
                            subgraph.add_node(n1, value=[n[1] for n in graph.nodes if n[0] == n1][0])
                            subgraph.add_edge(n1, n2, value)
    
    return subgraph

# Create the graph
graph = Graph()

graph.add_node(1, "C")
graph.add_node(2, "C")
graph.add_node(3, "C")
graph.add_node(4, "C")
graph.add_node(5, "C")
graph.add_node(6, "C")
graph.add_node(7, "C")

name = ["meth", "eth", "prop", "but", "pent", "hex", "hept", "oct", "non", "dec",\
            "undec", "dodec", "tridec", "tetradec", "pentadec", "hexadec", "heptadec",\
            "octadec", "nonadec", "icos"]

name_2 = ["di", "tri"]

def process(graph, depth=0, suffix=0, start_node=None):
    global name
    global bond_info
    global buffer
    global name_2
    if graph.nodes == []:
        return
    
    all_chain_output = all_chain_all(graph, start_node)
    
    def double_bond(graph, chain):
        double = []
        triple = []
        i = 0
        while i < len(chain[0]):
            for edge in graph.edges:
                if chain[0][i] in edge[0]:
                    if edge[1] == 2:
                        double.append(i+1)
                        i = i + 1
                        break
                    elif edge[1] == 3:
                        triple.append(i+1)
                        i = i + 1
                        break
            i = i + 1
        return chain + [double, triple]
    
    all_chain_output = [double_bond(graph, x) for x in all_chain_output]
    
    max_double = max([len(x[2])+len(x[3]) for x in all_chain_output])
    if max_double > 0:
        all_chain_output = [x for x in all_chain_output if max_double == len(x[2])+len(x[3])]
        least_sum = min([sum(x[2])+sum(x[3]) for x in all_chain_output])
        all_chain_output = [x for x in all_chain_output if least_sum == sum(x[2])+sum(x[3])]
    max_len = max([len(x[0]) for x in all_chain_output])
    
    all_chain_output = [x for x in all_chain_output if max_len == len(x[0])]
    all_chain_output = sorted(all_chain_output, key=lambda x: sum(x[1]))[0]

    sub_chain = []
    
    for i in list(set(all_chain_output[1])):
        count = all_chain_output[1].count(i)
        sel_edge = []
        for item in [x[0] for x in graph.edges]:
            if all_chain_output[0][i-1] == item[0]:
                sel_edge.append(item[1])
            elif all_chain_output[0][i-1] == item[1]:
                sel_edge.append(item[0])
        sel_edge = list(set(sel_edge) - set(all_chain_output[0]))
        
        for j in range(count):
            forbidden = list(set(all_chain_output[0] + sel_edge) - set([sel_edge[j]]))
                
            sub_chain.append([str(i), process(explorable_subgraph(graph, sel_edge[j], forbidden), depth+1, i, sel_edge[j])])
    output = []
    sub_chain_type = {}
    for item in sub_chain:
        if item[1] in sub_chain_type.keys():
            sub_chain_type[item[1]].append(item[0])
        else:
            sub_chain_type[item[1]] = [item[0]]
    for key in sorted(sub_chain_type.keys()):
        if len(sub_chain_type[key]) > 1:
            output.append(",".join(sub_chain_type[key]) + "-" + name_2[len(sub_chain_type[key])-2] + "(" + key + ")")
        else:
            output.append(",".join(sub_chain_type[key]) + "-(" + key + ")")
    name_list = [name[max_len-1]]
    if depth == 0:
        is_bond = False
        for t, item in enumerate(all_chain_output[2:4]):
            if len(item) == 1:
                name_list.append("-" + ",".join([str(x) for x in item]) + "-" + ["ene","yne"][t])
                is_bond = True
            elif len(item) > 1:
                name_list.append("-" + ",".join([str(x) for x in item]) + "-" + name_2[len(item)-2] + ["ene","yne"][t])
                is_bond = True
        if is_bond is False:
            name_list.append("ane")
    else:
        is_bond = False
        if is_bond is False:
            name_list.append("yl")
    return "-".join(output) + "".join(name_list)

def is_connected(edge_list):
    if not edge_list:
        return True

    graph = defaultdict(list)
    nodes = set()
    
    for u, v in edge_list:
        graph[u].append(v)
        graph[v].append(u)
        nodes.add(u)
        nodes.add(v)
    
    start_node = next(iter(nodes))
    visited = set()
    
    def bfs(node):
        queue = deque([node])
        visited.add(node)
        
        while queue:
            current = queue.popleft()
            for neighbor in graph[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
    
    bfs(start_node)
    
    return len(visited) == len(nodes)

def has_cycle(graph, node, visited, parent):
    visited.add(node)
    for neighbor in graph[node]:
        if neighbor not in visited:
            if has_cycle(graph, neighbor, visited, node):
                return True
        elif parent is not None and neighbor != parent:
            return True
    return False

def is_tree_with_max_4_children(edge_list):
    if not edge_list:
        return False  # An empty edge list cannot represent a tree.
    
    graph = defaultdict(list)
    nodes = set()
    
    for u, v in edge_list:
        graph[u].append(v)
        graph[v].append(u)
        nodes.add(u)
        nodes.add(v)
    if len(nodes) != 7:
        return False
    # Check connectivity
    if not is_connected(edge_list):
        return False
    
    # Check for cycles
    visited = set()
    start_node = next(iter(nodes))
    if has_cycle(graph, start_node, visited, None):
        return False
    
    # Check the number of edges
    if len(edge_list) != len(nodes) - 1:
        return False
    
    # Check for maximum 4 children per node
    for node in graph:
        if len(graph[node]) > 4:
            return False
        
    return True
def generate_tree(num):
    node_pair = [list(x) for x in itertools.combinations(range(1,num+1),2)]
    output = []
    for i in range(1,len(node_pair)):
        for item in itertools.combinations(node_pair, i):
            item = [list(x) for x in item]
            if is_tree_with_max_4_children(item):
                output.append(item)
    return output
final = []
for new in generate_tree(7):
    #new = generate_constrained_graph_5nodes()
    graph.edges = [[[y[0], y[1]], 1] for y in new]
    x = process(graph).replace("-1-", "")
    for item in name:
        x = x.replace("("+item+"yl)", item+"yl")
    final.append(x)
final = list(set(final))
for item in final:
    print(item)
add_hydrogens(graph)
#print(graph.nodes)
#print(graph.edges)
#print(graph)
# Draw the graph
draw_graph(graph)
