from collections import defaultdict
import matplotlib.pyplot as plt
from igraph import Graph as IGraph, plot
import itertools
import copy
class Graph:
    # Atom valency information for common elements (neutral charge)
    atom_valency = {
        "C": 4,  # Carbon forms 4 bonds
        "O": 6,  # Oxygen forms 2 bonds
        "N": 5,  # Nitrogen forms 3 bonds
        "H": 1   # Hydrogen forms 1 bond
    }
    def __init__(self):
        self.nodes = []  # List to store nodes
        self.edges = []  # List to store edges
    def add_node(self, atom_type="C", charge_lp=[0,0]):
        """Add a node with atom type and charge/lone pairs."""
        self.nodes.append([atom_type, charge_lp.copy()])  # Store atom type and charge/lone pairs

    def add_edge(self, node_1, node_2, bond=1):
        """Add an edge between two nodes with a specified bond (1 = single, 2 = double, etc.)."""
        if node_1 > node_2:
            node_1, node_2 = node_2, node_1
        self.edges.append([node_1, node_2, bond])

    def get_bond_type(self, node_1, node_2):
        """Get the bond type between two nodes."""
        if node_1 > node_2:
            node_1, node_2 = node_2, node_1
        for item in self.edges:
            if item[0] == node_1 and item[1] == node_2:
                return item[2]
        return None
    def canonical_label(self):
        """Generate a canonical label for the graph."""
        # Step 1: Sort nodes by their properties and remember the new indices
        sorted_nodes_with_indices = sorted(
            (i, (node[0], tuple(node[1]))) for i, node in enumerate(self.nodes)
        )  # Ensure 'charge_lp' is a tuple
        node_mapping = {old: new for new, (old, _) in enumerate(sorted_nodes_with_indices)}

        # Step 2: Remap edges based on the new indices
        remapped_edges = [
            (node_mapping[node_1], node_mapping[node_2], bond)
            for node_1, node_2, bond in self.edges
        ]

        # Step 3: Sort edges lexicographically
        sorted_edges = sorted(remapped_edges)

        # Step 4: Convert everything to tuples to make it hashable
        return (
            tuple(node[1] for node in sorted_nodes_with_indices),  # Canonical node list
            tuple(sorted_edges)  # Canonical edge list
        )


    def count_bonds(self, node_index):
        """Count the total number of bonds for a given node index."""
        bond_count = 0
        for edge in self.edges:
            if edge[0] == node_index or edge[1] == node_index:
                bond_count += edge[2]  # Increment bond count by bond type (1, 2, etc.)
        return bond_count

    def update_lone_pairs(self):
        """Update lone pairs for all nodes based on their current bonding situation."""
        for node_index in range(len(self.nodes)):
            atom_type, charge = self.nodes[node_index]
            if atom_type is None:
                continue
            current_bond_count = self.count_bonds(node_index)
            hydrogens_needed = 0
            
            if (self.atom_valency[atom_type] - charge[0]) >= 4:
                hydrogens_needed = 8 - ((self.atom_valency[atom_type] - charge[0]) + current_bond_count)
            else:
                hydrogens_needed = self.atom_valency[atom_type] - charge[0] - current_bond_count
            lone_pair_count = self.atom_valency[atom_type] - charge[0] - current_bond_count - hydrogens_needed
            self.nodes[node_index][1][1] = int(lone_pair_count/2)  # Update lone pair count

    def add_hydrogens(self):
        """Automatically add hydrogens to the molecule based on atoms' valency and current bonds."""
        orig = copy.deepcopy(self.nodes)
        for i, (atom_type, charge) in enumerate(orig):

            # Calculate current bond count
            current_bond_count = self.count_bonds(i)
            hydrogens_needed = 0
            if (self.atom_valency[atom_type] - charge[0]) >= 4:
                hydrogens_needed = 8 - ((self.atom_valency[atom_type] - charge[0]) + current_bond_count)
            else:
                hydrogens_needed = self.atom_valency[atom_type] - charge[0] - current_bond_count
            # Ensure we only add hydrogens if more are needed
            for _ in range(hydrogens_needed):
                hydrogen_index = len(self.nodes)  # The new node index for the hydrogen
                self.add_node("H")  # Add a hydrogen atom
                self.add_edge(i, hydrogen_index)  # Connect the hydrogen atom with a single bond


def apply(compound, formula, formula_rhs):
    output = []
    for item in itertools.permutations(enumerate(compound.nodes), len(formula.nodes)):
        if not all(formula.get_bond_type(x[0][0], x[1][0]) == compound.get_bond_type(x[0][1][0], x[1][1][0]) for x in itertools.combinations(enumerate(item), 2)):
            continue
        failed = False
        data = {}
        for i in range(len(item)):
            if formula.nodes[i][0][:2] == "u_":
                if formula.nodes[i][0] not in data.keys():
                    data[formula.nodes[i][0]] = item[i][0]
                else:
                    if data[formula.nodes[i][0]] != item[i][0]:
                        failed = True
                        break
            else:
                if formula.nodes[i][0] != item[i][1][0]:
                    failed = True
                    break
            if formula.nodes[i][1] != item[i][1][1]:
                failed = True
                break
        if failed == False:
            print(data)
            output.append(copy.deepcopy(compound))
            for i in range(len(item)):
                output[-1].nodes[item[i][0]] = copy.deepcopy(formula_rhs.nodes[i])
                if formula_rhs.nodes[i][0] in data.keys():
                    output[-1].nodes[item[i][0]][0] = compound.nodes[data[formula_rhs.nodes[i][0]]][0]
            for x in itertools.combinations(enumerate(item), 2):
                a, b = x[0][1][0], x[1][1][0]
                if a> b:
                    a, b = b, a
                for i in range(len(output[-1].edges)):
                    if output[-1].edges[i][0] == a and output[-1].edges[i][1] == b:
                        output[-1].edges[i][2] = formula_rhs.get_bond_type(x[0][0], x[1][0])
    return output

def display(compound, flag=True):
    g = copy.deepcopy(compound)
    if flag:
        g.add_hydrogens()
    ig = IGraph()
    ig.add_vertices(len(g.nodes))  # Add vertices based on number of nodes
    ig.add_edges([(edge[0], edge[1]) for edge in g.edges])  # Add edges

    # Vertex and edge properties
    ig.vs["label"] = [str(node[1]).replace(" ","")[1:-1] if node[1] != [0,0] else "" for node in g.nodes]  # Show charge only if non-zero
    ig.es["bond"] = [edge[2] for edge in g.edges]  # Bond types

    # Define colors for different atom types
    atom_colors = {
        "C": "gray",
        "O": "red",
        "N": "green",
        "H": "white",
        "F": "purple"
    }
    vertex_colors = [atom_colors[node[0]] for node in g.nodes]

    # Edge colors for bonds
    bond_colors = ["blue" if bond == 1 else "red" if bond == 2 else "green" for bond in ig.es["bond"]]

    # Layout and plot settings
    layout = ig.layout("auto")
    visual_style = {
        "vertex_label": ig.vs["label"],  # Show only the charge inside circles if non-zero
        "vertex_size": 40,
        "vertex_color": vertex_colors,  # Color based on the atom type
        "edge_width": [2 if bond == 2 else 3 if bond == 3 else 1 for bond in ig.es["bond"]],  # Adjust width for different bonds
        "edge_color": bond_colors,  # Color based on bond type
        "layout": layout,
        "bbox": (300, 300),
        "margin": 20
    }

    # Create a matplotlib figure
    fig, ax = plt.subplots(figsize=(6, 6))

    # Plot the igraph object
    plot(ig, target=ax, **visual_style)

    # Show the plot
    plt.show()

def are_chemicals_same(graph1, graph2):
    """Compare if two ChemicalGraph objects represent the same compound."""
    return graph1.canonical_label() == graph2.canonical_label()
def count_chemical_graphs(graph_list):
    """Count occurrences of unique chemical graphs."""
    counter = defaultdict(int)
    
    for graph in graph_list:
        # Get the canonical label for the graph
        canonical = graph.canonical_label()
        # Increment the count for this canonical label
        counter[canonical] += 1
    
    return counter
def find_edge(lst, node1, node2):
    for item in lst:
        if node1 in item[:2] and node2 in item[:2]:
            return item[2]
    return None
def check(var_list, a, b):
    orig = copy.deepcopy(var_list)
    for i in range(len(a)):
        if type(b[i])==str and "v_" in b[i]:
            if b[i] in var_list.keys():
                if var_list[b[i]] != a[i]:
                    return False, orig
            else:
                var_list[b[i]] = a[i]
        elif a[i] != b[i]:
            return False, orig
    return True, orig
def flat(lst):
    output = []
    for item in lst:
        output += [item[0], item[1][0], item[1][1]]
    return output

def applyformula(compound, formula, nxt):
    output = []
    for sub in itertools.permutations(enumerate(compound.nodes), len(formula.nodes)):
        if not all(find_edge(compound.edges, sub[i][0], sub[i+1][0]) is not None and find_edge(compound.edges, sub[i][0], sub[i+1][0])==find_edge(formula.edges, i, i+1) for i in range(len(sub)-1)):
            continue
        sub2 = [x[1] for x in sub]
        
        result = check({}, flat(sub2), flat(formula.nodes))
        if not result[0]:
            continue
        nxt2 = copy.deepcopy(nxt)
        compound2 = copy.deepcopy(compound)
        for key in result[1].keys():
            for i in range(len(nxt2.nodes)):
                if nxt2.nodes[i][0] is not None:
                    nxt2.nodes[i][0] = nxt2.nodes[i][0].replace(key, result[1][key])
                if nxt2.nodes[i][1] is not None:
                    nxt2.nodes[i][1] = nxt2.nodes[i][1].replace(key, result[1][key])
                if nxt2.nodes[i][1][0] is not None:
                    nxt2.nodes[i][1][0] = nxt2.nodes[i][1][0].replace(key, result[1][key])
                if nxt2.nodes[i][1] is not None:
                    nxt2.nodes[i][1][1] = nxt2.nodes[i][1][1].replace(key, result[1][key])
        
        for i in range(len(nxt2.nodes)):
            compound2.nodes[sub[i][0]] = copy.deepcopy(nxt2.nodes[i])
            
        for i in range(3):
            for x in nxt2.edges:
                if i == x[2]:
                    for j in range(len(compound2.edges)):
                        if sub[x[0]][0] in compound2.edges[j][:2] and sub[x[1]][0] in compound2.edges[j][:2]:
                            compound2.edges[j][2] = x[2]

        enter = True
        while enter:
            enter = False
            for i in range(len(compound2.nodes)-1,-1,-1):
                if compound2.nodes[i][0] is None:
                    for j in range(len(compound2.edges)-1,-1,-1):
                        if i in compound2.edges[j][:2]:
                            compound2.edges.pop(j)
                    for j in range(len(compound2.edges)-1,-1,-1):
                        for k in range(2):
                            if compound2.edges[j][k] > i:
                                compound2.edges[j][k] = compound2.edges[j][k] - 1
                    
                    
                    
                    compound2.nodes.pop(i)
                    
                    enter = True
                    break
        output.append(copy.deepcopy(compound2))
    return output
a = Graph()
a.add_node("C")
a.add_node("C", [1, 0])
a.add_edge(0, 1)
a.update_lone_pairs()
a.add_hydrogens()

b = Graph()
b.add_node("H")
b.add_node("C")
b.add_node("C", [1, 0])
b.add_edge(0, 1)
b.add_edge(1, 2)
b.update_lone_pairs()

c = Graph()
c.add_node(None)
c.add_node("C")
c.add_node("C")
c.add_edge(0, 1, 1)
c.add_edge(1, 2, 2)
c.update_lone_pairs()

a = Graph()
a.add_node("C")
a.add_node("C")
a.add_node("C")
a.add_edge(0, 1)
a.add_edge(1, 2, 2)
a.update_lone_pairs()
a.add_hydrogens()


b = Graph()
b.add_node("H")
b.add_node("C")
b.add_node("C")
b.add_node("C")
b.add_edge(0, 1)
b.add_edge(1, 2)
b.add_edge(2, 3, 2)
b.update_lone_pairs()

c = Graph()
c.add_node(None)
c.add_node("C")
c.add_node("C")
c.add_node("C", [-1, 0])
c.add_edge(0, 1)
c.add_edge(1, 2, 2)
c.add_edge(2, 3)
c.update_lone_pairs()

a = Graph()
a.add_node("C")
a.add_node("C")
a.add_node("C", [1, 0])
a.add_node("C")
a.add_node("C")
a.add_node("C")
a.add_edge(0, 1)
a.add_edge(1, 2)
a.add_edge(2, 3)
a.add_edge(3, 4)
a.add_edge(4, 0)
a.add_edge(2, 5)
a.update_lone_pairs()
a.add_hydrogens()

display(a)

b = Graph()
b.add_node("H")
b.add_node("C")
b.add_node("C", [1, 0])
b.add_edge(0, 1)
b.add_edge(1, 2)
b.update_lone_pairs()

c = Graph()
c.add_node(None)
c.add_node("C")
c.add_node("C")
c.add_edge(0, 1, 1)
c.add_edge(1, 2, 2)
c.update_lone_pairs()

result = applyformula(a, b, c)
result = count_chemical_graphs(result)
print(result)
eewwd()
def resonance(compound, depth=2):
    if depth == 0:
        return []
    compound.update_lone_pairs()
    output = [compound]
    
    a = Graph()
    a.add_node("O", [0,2])
    a.add_node("u_0")
    a.add_node("u_1", [-1,0])
    a.add_edge(0, 1, 2)
    a.add_edge(1, 2, 1)

    b = Graph()
    b.add_node("O", [-1,3])
    b.add_node("u_0")
    b.add_node("u_1")
    b.add_edge(0, 1, 1)
    b.add_edge(1, 2, 2)

    c = Graph()
    c.add_node("u_0", [1,-1])
    c.add_node("u_1", [0,1])
    c.add_edge(0, 1)

    d = Graph()
    d.add_node("u_0")
    d.add_node("u_1", [1,0])
    d.add_edge(0, 1, 2)
    
    for item in apply(compound, a, b) + apply(compound, b, a) + apply(compound, c, d) + apply(compound, d, c):
        output += resonance(item, depth-1)
    return output




g = Graph()
g.add_node()
g.add_node()
g.add_node()
g.add_node("C", [1,0])
g.add_node()
g.add_node("N")
g.add_node()
g.add_node()
g.add_edge(0, 1)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(3, 4)
g.add_edge(3, 5)
g.add_edge(5, 6)
g.add_edge(5, 7)

# Creating the graph
g = Graph()
g.add_node()
g.add_node()
g.add_node("C", [-1,0])
g.add_node("O")
g.add_edge(0, 1)
g.add_edge(1, 2)
g.add_edge(1, 3, 2)

for item in resonance(g):
    display(item)

