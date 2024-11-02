import itertools
import igraph as ig
import matplotlib.pyplot as plt
import copy
import random
from collections import defaultdict, deque
class Graph:
    def __init__(self):
        self.nodes = []  # Each node is a list [node_id, value]
        self.edges = []  # Each edge is a list [[node1_id, node2_id], value]
    def add_node(self, node_id, value="C"):
        if node_id not in [node[0] for node in self.nodes]:
            self.nodes.append([node_id, value])
    def edge_exists(self, node1_id, node2_id):
        return any(edge for edge, _ in self.edges if set(edge) == {node1_id, node2_id})
    def add_edge(self, node1_id, node2_id, value=1):
        self.edges.append([[node1_id, node2_id], value])

    def build_adjacency_list(self):
        adj_list = {node[0]: [] for node in self.nodes}
        for edge, value in self.edges:
            node1_id, node2_id = edge
            adj_list[node1_id].append(node2_id)
            adj_list[node2_id].append(node1_id)  # Ensure bidirectional edge
        return adj_list
    def degree(self, node_id):
        adj_list = self.build_adjacency_list()
        return len(adj_list.get(node_id, []))
    def find_all_cycles(self):
        adj_list = self.build_adjacency_list()
        visited = set()
        recStack = []
        all_cycles = []

        def dfs(v, parent):
            visited.add(v)
            recStack.append(v)

            for neighbor in adj_list[v]:
                if neighbor not in visited:
                    if dfs(neighbor, v):
                        return True
                elif neighbor != parent and neighbor in recStack:
                    # Found a cycle
                    cycle_start_index = recStack.index(neighbor)
                    cycle = recStack[cycle_start_index:] + [neighbor]
                    all_cycles.append(cycle)

            recStack.pop()
            return False

        for node_id, _ in self.nodes:
            if node_id not in visited:
                dfs(node_id, None)

        # Determine return value based on number of cycles found
        if len(all_cycles) == 0:
            return []  # No cycles found
        elif len(all_cycles) == 1:
            return list(set(all_cycles[0]))  # Exactly one cycle found
        else:
            return None  # More than one cycle found
    def get_edges(self, node_id):
        """Returns the edges connected to a specific node_id."""
        return [edge for edge in self.edges if node_id in edge[0]]
    def add_hydrogens(self):
        """Automatically add hydrogen atoms to carbon atoms to ensure each carbon has four bonds."""
        hydrogens = []  # To keep track of hydrogen nodes added
        count_hydrogen = len(self.nodes)+1
        #print(count_hydrogen)
        #print(self.find_all_cycles())
        for node_id, value in self.nodes:
            if value == "C":
                # Count existing bonds
                existing_bonds = sum([edge[1] for edge in self.get_edges(node_id)])
                needed_hydrogens = 4 - existing_bonds

                # Add hydrogens as needed
                for _ in range(needed_hydrogens):
                    hydrogen_id = count_hydrogen
                    count_hydrogen += 1
                    self.add_node(hydrogen_id, value="H")
                    self.add_edge(node_id, hydrogen_id)
                    hydrogens.append(hydrogen_id)

        return hydrogens
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
count_hydrogen = 0
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
    ig.plot(g, layout=layout, vertex_label=g.vs['label'], vertex_size=15, vertex_color="lightblue", edge_color=edge_colors, edge_width=2, bbox=(800, 600), margin=50, target="graph_igraph.png")
    
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

name = ["meth", "eth", "prop", "but", "pent", "hex", "hept", "oct", "non", "dec",\
            "undec", "dodec", "tridec", "tetradec", "pentadec", "hexadec", "heptadec",\
            "octadec", "nonadec", "icos"]
name =[
    "meth", "eth", "prop", "but", "pent", "hex", "hept", "oct", "non", "dec",
    "undec", "dodec", "tridec", "tetradec", "pentadec", "hexadec", "heptadec",
    "octadec", "nonadec", "icos", "henicos", "docos", "tricos", "tetracos",
    "pentacos", "hexacos", "heptacos", "octacos", "nonacos", "triacont",
    "tetracont", "pentacont", "hexacont", "heptacont", "octacont", "nonacont",
    "hect", "dodecacos", "tricosacont", "tetracontad", "pentacontad"
]

name_2 = [
    "di", "tri", "tetra", "penta", "hexa", "hepta", "octa", "nona", "deca", 
    "undeca", "dodeca", "trideca", "tetradeca", "pentadeca", "hexadeca", 
    "heptadeca", "octadeca", "nonadeca", "icosa", "henicosa", "docosa", 
    "tricosa", "tetracosa", "pentacosa", "hexacosa", "heptacosa", 
    "octacosa", "nonacosa", "triaconta", "tetraconta", "pentaconta", 
    "hexaconta", "heptaconta", "octaconta", "nonaconta"
]



def neighbor(graph, node):
    output = []
    for item in graph.edges:
        if node == item[0][0]:
            output.append(item[0][1])
        elif node == item[0][1]:
            output.append(item[0][0])
    return list(set(output))
            

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
            i = i + 1
        i = 0
        while i < len(chain[0]):
            for edge in graph.edges:
                if chain[0][i] in edge[0]:
                    if edge[1] == 3:
                        triple.append(i+1)
                        i = i + 1
                        break
            i = i + 1
        return chain + [double] + [triple]
    
    all_chain_output = [double_bond(graph, x) for x in all_chain_output]
    
    
    for i in range(len(all_chain_output)):
        all_chain_output[i] += [0]

    all_chain_new = []
    
    if graph.find_all_cycles() != []:
        cycle = graph.find_all_cycles()
        #print(cycle)
        for i in range(len(cycle)):
            new_cycle = [cycle[(i+j)%len(cycle)] for j in range(len(cycle))]
            #print(new_cycle)
            subchain = []
            for index, node in enumerate(new_cycle):
                subchain += [index+1]*(len(neighbor(graph, node))-2)
            #print(subchain)
            all_chain_output.append(double_bond(graph, [new_cycle, subchain])+[1])
                
    #for item in all_chain_output:
    #    print(item)
    
    max_double = max([len(x[2])+len(x[3]) for x in all_chain_output])
    if max_double > 0:
        all_chain_output = [x for x in all_chain_output if max_double == len(x[2])+len(x[3])]
        least_sum = min([sum(x[2])+sum(x[3]) for x in all_chain_output])
        all_chain_output = [x for x in all_chain_output if least_sum == sum(x[2])+sum(x[3])]

    max_cycle = max([x[4] for x in all_chain_output])
    all_chain_output = [x for x in all_chain_output if max_cycle == x[4]]
    
    max_len = max([len(x[0]) for x in all_chain_output])
    
    all_chain_output = [x for x in all_chain_output if max_len == len(x[0])]
    all_chain_output = sorted(all_chain_output, key=lambda x: sum(x[1]))
    all_chain_output = sorted(all_chain_output, key=lambda x: sum(x[2]))[0]

    sub_chain = []
    #print(all_chain_output)
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
    if all_chain_output[4] == 1:
        name_list = ["cyclo"]+name_list
    if depth == 0:
        if len(all_chain_output[2])!=0 and len(all_chain_output[3])!=0:
            arr = ["en", "yne"]
        else:
            arr = ["en", "yne"]
        is_bond = False
        for t, item in enumerate(all_chain_output[2:4]):
            if len(item) == 1:
                name_list.append("-" + ",".join([str(x) for x in item]) + "-" + arr[t])
                is_bond = True
            elif len(item) > 1:
                name_list.append("-" + ",".join([str(x) for x in item]) + "-" + name_2[len(item)-2] + arr[t])
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
        return True  # If there are no edges, consider the graph connected only if there are no nodes or one node.

    graph = defaultdict(list)
    nodes = set()

    # Build the graph
    for u, v in edge_list:
        graph[u].append(v)
        graph[v].append(u)
        nodes.add(u)
        nodes.add(v)

    # Pick an arbitrary start node
    start_node = next(iter(nodes))
    visited = set()

    # BFS to check connectivity
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

    # If we visited all nodes, the graph is connected
    return len(visited) == len(nodes)

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
    
    # Check connectivity
    if not is_connected(edge_list):
        return False
    
    # Check the number of edges (N-1 edges for a tree with N nodes)
    if len(nodes) != 6:
        return False
    
    # Check for maximum 4 children per node
    for node in graph:
        if len(graph[node]) > 4:
            return False
def generate_random_graph(n):
    graph = Graph()
    
    # Add nodes to the graph
    for i in range(n):
        graph.add_node(i)
    
    # Prim's algorithm to generate a random spanning tree
    connected_nodes = set()
    
    # Start with the first node
    current_node = random.randint(0, n - 1)
    connected_nodes.add(current_node)
    
    while len(connected_nodes) < n:
        possible_edges = []
        for node in connected_nodes:
            if graph.degree(node) < 4:  # Node has less than 4 neighbors
                for neighbor in range(n):
                    if (neighbor not in connected_nodes and
                            not graph.edge_exists(node, neighbor) and
                            graph.degree(neighbor) < 4):
                        possible_edges.append((node, neighbor))
        
        if possible_edges:
            edge = random.choice(possible_edges)
            graph.add_edge(edge[0], edge[1])
            connected_nodes.add(edge[1])
    
    # Optionally add one edge to create exactly one cycle, if needed
    potential_edges = list(itertools.combinations(range(n), 2))
    random.shuffle(potential_edges)
    
    for edge in potential_edges:
        if (not graph.edge_exists(edge[0], edge[1]) and
                graph.degree(edge[0]) < 4 and
                graph.degree(edge[1]) < 4):
            graph.add_edge(edge[0], edge[1])
            # Stop after adding one cycle
            break
    
    return graph
def generate_trees(num):
    if num <= 0:
        return []
    """
    node_pair = [list(x) for x in itertools.combinations(range(1, num + 1), 2)]
    output = []
    
    for i in range(1, len(node_pair)):
        for item in itertools.combinations(node_pair, i):
            item = [list(x) for x in item]
            if is_tree_with_max_4_children(item) and edge_convert(num, item).find_all_cycles() is not None:
                output.append(item)
    """
    return [generate_random_graph(40)]

    return item
def edge_convert(num, lst):
    graph = Graph()
    for i in range(num):
        graph.add_node(i+1)
    for i in range(len(lst)):
        graph.add_edge(lst[i][0], lst[i][1])
    return graph

# Basic data structure, which can nest to represent math equations
class TreeNode:
    def __init__(self, name, children=None):
        self.name = name
        self.children = children or []

# convert string representation into tree
def tree_form(tabbed_strings):
    lines = tabbed_strings.split("\n")
    root = TreeNode("Root") # add a dummy node
    current_level_nodes = {0: root}
    stack = [root]
    for line in lines:
        level = line.count(' ') # count the spaces, which is crucial information in a string representation
        node_name = line.strip() # remove spaces, when putting it in the tree form
        node = TreeNode(node_name)
        while len(stack) > level + 1:
            stack.pop()
        parent_node = stack[-1]
        parent_node.children.append(node)
        current_level_nodes[level] = node
        stack.append(node)
    return root.children[0] # remove dummy node

# convert tree into string representation
def str_form(node):
    def recursive_str(node, depth=0):
        result = "{}{}".format(' ' * depth, node.name) # spacings
        for child in node.children:
            result += "\n" + recursive_str(child, depth + 1) # one node in one line
        return result
    return recursive_str(node)

def replace(equation, find, r):
    if str_form(equation) == str_form(find):
      return r
    col = TreeNode(equation.name, [])
    for child in equation.children:
      col.children.append(replace(child, find, r))
    return col

def remove_past(equation):
    coll = TreeNode(equation.name, [])
    for child in equation.children:
      if child.name == "del":
          for subchild in child.children:
              coll.children.append(remove_past(subchild))
      else:
          coll.children.append(remove_past(child))
    return coll

def break_equation(equation):
    sub_equation_list = [equation]
    equation = equation
    for child in equation.children: # breaking equation by accessing children
        sub_equation_list += break_equation(child) # collect broken equations
    return sub_equation_list

def multiple(equation):
    def split(equation):
        output = []
        for i in range(len(equation.children)-1):
            output.append(TreeNode("compound", [equation.children[i], equation.children[-1]]))
        return TreeNode("del", output)
    for item in break_equation(equation):
        if item.name == "compound" and len(item.children) > 2:
            equation = replace(equation, item, split(item))
    return remove_past(equation)
    #return equation

def pre(compound):
    global name_2
    
    compound = compound.replace("ane", "")

    for item in name:
        if item == compound:
            return "0*"+compound
        elif "cyclo"+item == compound:
            return "0*"+compound
    
    for i in range(10,-1,-1):
        compound = compound.replace(str(i)+"-", str(i)+"*")
        
    for item in name_2:
        compound = compound.replace(item, "")
        
    compound = compound.replace("yl", ";")
    
    #compound = compound.replace("cyclo", "cyclo+")
    compound = compound.replace(";)", ")")
    compound = compound.replace(";;", ";")
    compound = compound.replace(")cyclo", ");cyclo")
    compound = compound.replace(";cyclo", ";cyclo")

    for item in name:
        compound = compound.replace(";" + item, ";0*" + item)
    compound = compound.replace(";-", ";")
    compound = compound.replace(";cyclo", ";0*cyclo")
    compound = compound.replace(")-", ");")
    return compound
def name2compound(graph, name_eq, position):
    global name
    if name_eq.name != "segment":
        return name2compound(copy.deepcopy(graph), TreeNode("segment", [name_eq]), position)
    for child in name_eq.children:
        if int(child.children[0].name) == 0:
            graph = add_base(copy.deepcopy(graph), name.index(child.children[1].name.replace("cyclo", ""))+1, position, child.children[1].name.replace("cyclo", "") != child.children[1].name)
    #draw_graph(graph)
    for child in name_eq.children:
        if int(child.children[0].name) == 0:
            continue
        if child.children[1].name != "segment":
            graph = add_base(copy.deepcopy(graph), name.index(child.children[1].name.replace("cyclo", ""))+1, position+int(child.children[0].name), child.children[1].name.replace("cyclo", "") != child.children[1].name)
        else:
            graph2 = name2compound(Graph(), child.children[1], 0)
            max_node = sorted(graph.nodes, key=lambda x: -x[0])[0][0]
            for i in range(len(graph2.nodes)):
                graph2.nodes[i][0] += max_node
            for i in range(len(graph2.edges)):
                graph2.edges[i][0][0] += max_node
                graph2.edges[i][0][1] += max_node
            graph.nodes += graph2.nodes
            graph.edges += graph2.edges
            graph.add_edge(position+int(child.children[0].name), max_node+1)
    return graph
"""
x = process(graph)
for item in name:
    x = x.replace("("+item+"yl)", item+"yl")
print(x)
"""
def add_base(graph, num, position, cyclic=False):
    max_node = 0
    if position != 0:
        max_node = sorted(graph.nodes, key=lambda x: -x[0])[0][0]
    
    for i in range(num):
        graph.add_node(max_node+i+1)
        if i != 0:
            graph.add_edge(max_node+i, max_node+i+1)
    if cyclic:
        graph.add_edge(max_node+1, max_node+num)
    if position != 0:
        graph.add_edge(position, max_node+1)
    return graph
import chemistry_parser_2
def post(x):
    global name
    for item in name:
        x = x.replace("("+item+"yl)", item+"yl")
    return x
final = []

# fancy print
def string_equation_helper(equation_tree):
    if equation_tree.children == []:
        return equation_tree.name # leaf node
    s = "(" # bracket
    if len(equation_tree.children) == 1:
        s = equation_tree.name + s
    #sign = {"f_add": "+", "f_mul": "*", "f_pow": "^", "f_div": "/", "f_int": ",", "f_sub": "-", "f_dif": "?", "f_sin": "?", "f_cos": "?", "f_tan": "?", "f_eq": "=", "f_sqt": "?"} # operation symbols
    sign = {"compound": "-", "segment": "?"}
    for child in equation_tree.children:
        s+= string_equation_helper(copy.deepcopy(child)) + sign[equation_tree.name]
    s = s[:-1] + ")"
    return s

# fancy print main function
def string_equation(eq): 
    eq = eq.replace("v_0", "x")
    eq = eq.replace("v_1", "y")
    eq = eq.replace("v_2", "z")
    eq = eq.replace("d_", "")
    
    return string_equation_helper(tree_form(eq))

x = "6-(2-methylpropyl)-7-ethyl-3,5-dimethyldodecane"
eq = chemistry_parser_2.take_input(pre(x))
eq = multiple(tree_form(eq))
print(string_equation(str_form(eq)))
g = name2compound(Graph(), eq, 0)

print(g.nodes, g.edges)
#print()
#draw_graph(g)
y = post(process(g))
print(x==y)

for new in generate_trees(6):
    #new = generate_constrained_graph_5nodes()
    #print(new)
    #tmp = edge_convert(6, new)
    #if tmp.find_all_cycles() != []:
    #    draw_graph(tmp)
    tmp = new
    #sdraw_graph(tmp)
    
    #print(new.find_all_cycles())
    x = process(tmp)
    
    for item in name:
        x = x.replace("("+item+"yl)", item+"yl")
    print(x)
    if x not in final:
        #tmp6 = x
        #x = "1,3,5-trimethylhexane"
        eq = chemistry_parser_2.take_input(pre(x))
        #print(eq)
        eq = multiple(tree_form(eq))
        print(string_equation(str_form(eq)))
        g = name2compound(Graph(), eq, 0)
        print(len(g.nodes))
        print(g.nodes)
        print(g.edges)
        g.add_hydrogens()
        draw_graph(g)
        count =0
        for item in g.nodes:
            if item[1] == "H":
                count += 1
        print(count)
        break
        #print()
        #draw_graph(g)
        y = post(process(g))
        #print(x==y)
        print(x)
        final.append(x)
    
#
#print(graph.nodes)
#print(graph.edges)
#print(graph)
# Draw the graph
#draw_graph(graph)

