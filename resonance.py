import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
import functools
from collections import defaultdict
import matplotlib.pyplot as plt
from igraph import Graph as IGraph, plot
import itertools
import copy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import numpy as np
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
    def add_node(self, atom_type="C", charge_lp=[0,0,0]):
        """Add a node with atom type and charge/lone pairs."""
        if len(charge_lp)==2:
            charge_lp += [0]
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
        def lists_to_tuples(nested_list):
            if isinstance(nested_list, list):
                return tuple(lists_to_tuples(item) for item in nested_list)
            return nested_list
        return lists_to_tuples([self.nodes, self.edges])


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
            if atom_type != "H":
                for edge in self.edges:
                    if node_index in edge[:2] and (self.nodes[edge[0]][0] == "H" or self.nodes[edge[1]][0] == "H"):
                        current_bond_count = current_bond_count - 1
            hydrogens_needed = 0
            
            if (self.atom_valency[atom_type] - charge[0]) >= 4:
                hydrogens_needed = 8 - ((self.atom_valency[atom_type] - charge[0]) + current_bond_count)
            else:
                hydrogens_needed = self.atom_valency[atom_type] - charge[0] - current_bond_count
            lone_pair_count = self.atom_valency[atom_type] - charge[0] - current_bond_count - hydrogens_needed - charge[2]
            self.nodes[node_index][1][1] = int(lone_pair_count/2)  # Update lone pair count
    def return_h_count(self):
        output = []
        for node_index in range(len(self.nodes)):
            atom_type, charge = self.nodes[node_index]
            if atom_type is None:
                continue
            current_bond_count = self.count_bonds(node_index)
            if atom_type != "H":
                for edge in self.edges:
                    if node_index in edge[:2] and (self.nodes[edge[0]][0] == "H" or self.nodes[edge[1]][0] == "H"):
                        current_bond_count = current_bond_count - 1
            hydrogens_needed = 0
            
            if (self.atom_valency[atom_type] - charge[0]) >= 4:
                hydrogens_needed = 8 - ((self.atom_valency[atom_type] - charge[0]) + current_bond_count)
            else:
                hydrogens_needed = self.atom_valency[atom_type] - charge[0] - current_bond_count
            output.append(hydrogens_needed)
        return tuple(output)
    def is_neutral(self):
        for _, charge in self.nodes:
            if charge[0] != 0:
                return False
        return True
    def pi_count(self):
        total = 0
        for _a, _b, x in self.edges:
            if x > 1:
                total += x - 1
        return total
    def add_hydrogens(self):
        """Automatically add hydrogens to the molecule based on atoms' valency and current bonds."""
        orig = copy.deepcopy(self)
        for i, (atom_type, charge) in enumerate(orig.nodes):

            # Calculate current bond count
            current_bond_count = orig.count_bonds(i)
            hydrogens_needed = 0
            if (self.atom_valency[atom_type] - charge[0]) >= 4:
                hydrogens_needed = 8 - ((self.atom_valency[atom_type] - charge[0]) + current_bond_count)
            else:
                hydrogens_needed = self.atom_valency[atom_type] - charge[0] - current_bond_count
            hydrogens_needed = hydrogens_needed - charge[2]
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

import copy
from io import BytesIO

def display(compound, flag=True):
    g = copy.deepcopy(compound)
    if flag:
        g.add_hydrogens()
    ig = IGraph()
    ig.add_vertices(len(g.nodes))  # Add vertices based on number of nodes
    ig.add_edges([(edge[0], edge[1]) for edge in g.edges])  # Add edges

    # Vertex and edge properties
    ig.vs["label"] = [
        str(node[1]).replace(" ", "")[1:-1] if node[1] != [0, 0, 0] else ""
        for node in g.nodes
    ]  # Show charge only if non-zero

    for i in range(len(ig.vs["label"])):
        if ig.vs[i]["label"][-2:] == ",0":
            ig.vs[i]["label"] = ig.vs[i]["label"][:-2]
    ig.es["bond"] = [edge[2] for edge in g.edges]  # Bond types

    # Define colors for different atom types
    atom_colors = {
        "C": "gray",
        "O": "red",
        "N": "green",
        "H": "white",
        "F": "purple",
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
        "margin": 20,
    }

    # Create a matplotlib figure
    fig, ax = plt.subplots(figsize=(6, 6))

    # Plot the igraph object
    plot(ig, target=ax, **visual_style)

    buf = BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)

    # Close the plot
    plt.close(fig)

    # Open the plot as a PIL Image
    img = Image.open(buf)

    return img


from PIL import Image
import math
def create_square_album(image_list, image_size=(300, 300)):
    # Resize images to the desired size
    images = [img.resize(image_size) for img in image_list]
    
    # Determine the grid dimensions for the square album
    num_images = len(images)
    grid_size = math.ceil(math.sqrt(num_images))  # Number of rows/columns in the square grid

    # Calculate the size of the album canvas
    album_size = (grid_size * image_size[0], grid_size * image_size[1])
    album = Image.new("RGB", album_size, color="white")
    
    # Paste images onto the canvas in a grid layout
    for idx, img in enumerate(images):
        row = idx // grid_size
        col = idx % grid_size
        x_offset = col * image_size[0]
        y_offset = row * image_size[1]
        album.paste(img, (x_offset, y_offset))

    return album



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

def check(var_list, a, b):
    orig = copy.deepcopy(var_list)
    for i in range(len(a)):
        if type(b[i])==str and b[i]=="nonzero":
            if a[i] == 0:
                return False, orig
        elif type(b[i])==str and b[i]=="any":
            continue
        elif type(b[i])==str and "v_" in b[i]:
            if b[i] in var_list.keys():
                
                if var_list[b[i]] != a[i]:
                    return False, orig
            else:
                
                var_list[b[i]] = a[i]
        elif a[i] != b[i]:
            return False, orig
        
    return True, var_list
def flat(lst):
    output = []
    for item in lst:
        if isinstance(item, list):  # Check if the item is a list
            output.extend(flat(item))  # Recursively flatten the sublist
        else:
            output.append(item)  # If it's not a list, just append the item
    return output
def convert_to_networkx(graph):
    G = nx.Graph()

    # Add nodes with their attributes
    for i, node in enumerate(graph.nodes):
        atom_type, charge = node
        G.add_node(i, atom_type=atom_type, charge=tuple(charge))

    # Add edges with bond attributes
    for edge in graph.edges:
        node1, node2, bond = edge
        G.add_edge(node1, node2, bond=bond)
    
    return G

def find_edge(lst, node1, node2):
    for item in lst:
        if node1 in item[:2] and node2 in item[:2]:
            return item[2]
    return None
def check2(var_list, a, b):
    for i in range(len(a)):
        if a[i] != b[i]:
            return False
    return True
def equalchem2(compound, formula):
    for sub in itertools.permutations(enumerate(compound.nodes), len(formula.nodes)):
        if not all(find_edge(compound.edges, sub[i][0], sub[j][0])==find_edge(formula.edges, i, j) for i,j in itertools.combinations(range(len(formula.nodes)), 2)):
            
            continue
        
        sub2 = [x[1] for x in sub]
        result = check2({}, flat(sub2), flat(formula.nodes))
        if not result[0]:
            continue
        return True
    return False
def equalchem(a, b):
    G1 = convert_to_networkx(a)
    G2 = convert_to_networkx(b)
    matcher = GraphMatcher(G1, G2, node_match=lambda n1, n2: n1['atom_type'] == n2['atom_type'] and n1['charge'] == n2['charge'], 
                       edge_match=lambda e1, e2: e1['bond'] == e2['bond'])
    return matcher.is_isomorphic()

def applyformula(compound, formula, nxt):
    output = []
    for sub in itertools.permutations(enumerate(compound.nodes), len(formula.nodes)):
        if not all(find_edge(compound.edges, sub[i][0], sub[j][0])==find_edge(formula.edges, i, j) for i,j in itertools.combinations(range(len(formula.nodes)), 2)):
            
            continue
        
        sub2 = [x[1] for x in sub]
        result = check({}, flat(sub2), flat(formula.nodes))
        if not result[0]:
            continue
        nxt2 = copy.deepcopy(nxt)
        compound2 = copy.deepcopy(compound)
        for key in result[1].keys():
            for i in range(len(nxt2.nodes)):
                if nxt2.nodes[i][0]==key:
                    nxt2.nodes[i][0] = result[1][key]
                if nxt2.nodes[i][1][0]==key:
                    nxt2.nodes[i][1][0] = result[1][key]
                if nxt2.nodes[i][1][1]==key:
                    nxt2.nodes[i][1][1] = result[1][key]
                if nxt2.nodes[i][1][2]==key:
                    nxt2.nodes[i][1][2] = result[1][key]
        
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
print()
done2 = []
done3 = []
def resonance(compound, depth=12):
    global done2
    global done3
    if depth == 0:
        return
    compound.update_lone_pairs()
    if compound.canonical_label() not in done3:
        done3.append(compound.canonical_label())
        
        done2.append(compound)
    else:
        return
    
    formula_a = []
    formula_b = []
    
    a = Graph()
    a.add_node("v_0",[0,"any"])
    a.add_node("v_1",[0,"any"])
    a.add_edge(0, 1, 2)
    
    b = Graph()
    b.add_node("v_0",[-1,0])
    b.add_node("v_1",[1,0])
    b.add_edge(0, 1)
    
    formula_a.append(a)
    formula_b.append(b)

    b = Graph()
    b.add_node("v_0",[-1,"any"])
    b.add_node("v_1",[1,"any"])
    b.add_edge(0, 1)
    
    a = Graph()
    a.add_node("v_0",[0,0])
    a.add_node("v_1",[0,0])
    a.add_edge(0, 1, 2)
    
    formula_a.append(b)
    formula_b.append(a)
    
    a = Graph()
    a.add_node("v_0",[1,"any"])
    a.add_node("v_1",[0,"any"])
    a.add_edge(0, 1)
    
    b = Graph()
    b.add_node("v_0",[0,0])
    b.add_node("v_1",[1,0])
    b.add_edge(0, 1, 2)
    
    formula_a.append(a)
    formula_b.append(b)

    b = Graph()
    b.add_node("v_0",[0,"any"])
    b.add_node("v_1",[1,"any"])
    b.add_edge(0, 1, 2)
    
    a = Graph()
    a.add_node("v_0",[1,0])
    a.add_node("v_1",[0,0])
    a.add_edge(0, 1)
    
    formula_a.append(b)
    formula_b.append(a)
    
    d = []
    for i in range(len(formula_a)):
        d += applyformula(compound, formula_a[i], formula_b[i])
    for item in d:
        if item.return_h_count() != compound.return_h_count():
            continue
        resonance(item, depth-1)

# Creating the graph


g = Graph()
g.add_node("C")
g.add_node("C")
g.add_node("C", [1,0])
g.add_edge(0, 1, 2)
g.add_edge(1, 2)
g.update_lone_pairs()
g.add_hydrogens()

g = Graph()
g.add_node()
g.add_node()
g.add_node("C", [-1,0])
g.add_node("O")
g.add_edge(0, 1)
g.add_edge(1, 2)
g.add_edge(1, 3, 2)
g.update_lone_pairs()
g.add_hydrogens()


g = Graph()
g.add_node("C")
g.add_node("C")
g.add_node("N")
g.add_edge(0, 1, 2)
g.add_edge(1, 2)
g.update_lone_pairs()
g.add_hydrogens()

g = Graph()
g.add_node("C")
g.add_node("C")
g.add_node("C", [-1,0])
g.add_edge(0, 1, 2)
g.add_edge(1, 2)
g.update_lone_pairs()
g.add_hydrogens()

g = Graph()
g.add_node()
g.add_node()
g.add_node()
g.add_node()
g.add_node("O")
g.add_edge(0, 1)
g.add_edge(1, 2)
g.add_edge(2, 3, 2)
g.add_edge(1, 4, 2)
g.update_lone_pairs()
g.add_hydrogens()

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
g.update_lone_pairs()
g.add_hydrogens()

g = Graph()
g.add_node()
g.add_node("N", [1,0])
g.add_node("O")
g.add_node("O",[-1,0])
g.add_edge(0, 1)
g.add_edge(1, 2, 2)
g.add_edge(1, 3)
g.update_lone_pairs()
g.add_hydrogens()

g = Graph()
g.add_node()
g.add_node()
g.add_node()
g.add_node()
g.add_node()
g.add_node()
g.add_node()
g.add_node("O")
g.add_edge(0, 1, 2)
g.add_edge(1, 2)
g.add_edge(2, 3, 2)
g.add_edge(3, 4)
g.add_edge(3, 5, 2)
g.add_edge(5, 0)
g.add_edge(0, 6)
g.add_edge(6, 7, 2)
g.update_lone_pairs()
g.add_hydrogens()

g = Graph()
g.add_node()
g.add_node()
g.add_node("O")
g.add_node()
g.add_node("O")
g.add_edge(0, 1)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(1, 4, 2)
g.update_lone_pairs()
g.add_hydrogens()

g = Graph()
g.add_node()
g.add_node()
g.add_node()
g.add_node()
g.add_node()
g.add_node()
g.add_node("O")
g.add_edge(0, 1, 2)
g.add_edge(1, 2)
g.add_edge(2, 3, 2)
g.add_edge(3, 4)
g.add_edge(4, 5, 2)
g.add_edge(5, 0)
g.add_edge(0, 6)
g.update_lone_pairs()
g.add_hydrogens()

g = Graph()
g.add_node()
g.add_node()
g.add_node()
g.add_node()
g.add_node()
g.add_node()
g.add_edge(0, 1, 2)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(3, 4)
g.add_edge(4, 5, 2)
g.add_edge(5, 0)
g.update_lone_pairs()
g.add_hydrogens()

g = Graph()
g.add_node()
g.add_node()
g.add_node("O")
g.add_edge(0, 1)
g.add_edge(1, 2, 2)
g.update_lone_pairs()
g.add_hydrogens()

g = Graph()
g.add_node("C")
g.add_node("C")
g.add_node("C")
g.add_node("C")
g.add_node("C")
g.add_node("O", [1,0])
g.add_edge(0, 1)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(3, 4)
g.add_edge(4, 5)
g.add_edge(5, 0, 2)
g.update_lone_pairs()
g.add_hydrogens()

g = Graph()
g.add_node()
g.add_node()
g.add_node("O")
g.add_node()
g.add_node("O")
g.add_edge(0, 1)
g.add_edge(1, 2)
g.add_edge(2, 3)
g.add_edge(1, 4, 2)
g.update_lone_pairs()
g.add_hydrogens()

g = Graph()
g.add_node("C",[-1,0,0])
g.add_node()
g.add_node()
g.add_node()
g.add_node("O", [1,0,0])
g.add_node()
g.add_edge(0, 1)
g.add_edge(1, 2, 2)
g.add_edge(2, 3)
g.add_edge(3, 4, 2)
g.add_edge(4, 5)
g.update_lone_pairs()
g.add_hydrogens()

image = []
resonance(g)
def cit(a, b):
    if a.is_neutral() and not b.is_neutral():
        return -1
    if not a.is_neutral() and b.is_neutral():
        return 1
    if a.pi_count() > b.pi_count():
        return -1
    if a.pi_count() < b.pi_count():
        return 1
    if len(a.nodes)==len(b.nodes):
        for i,j in itertools.combinations(range(len(a.nodes)), 2):
            origa = copy.deepcopy(a)
            origb = copy.deepcopy(b)
            origa.nodes[i] = ["C", [0,0,0]]
            origb.nodes[j] = ["C", [0,0,0]]
            if equalchem(origa, origb):
                if a.nodes[i][1][0] == b.nodes[j][1][0] and a.nodes[i][1][0] != 0:
                    e = {"C":6, "N":7, "O":8}
                    if a.nodes[i][0]*e[a.nodes[i][0]] > a.nodes[i][0]*e[b.nodes[j][0]]:
                        return -1
                    else:
                        return 1
        for i,j,m,n in itertools.permutations(range(len(a.nodes)), 4):
            if not(flat(a.nodes[i])==flat(b.nodes[m]) and flat(a.nodes[j])==flat(b.nodes[n]) and a.nodes[i][1][0] == 1 and\
               (a.nodes[j][1][1] > 0 or a.nodes[j][1][0] == -1)):
                continue
            origa = copy.deepcopy(a)
            origb = copy.deepcopy(b)
            origa.nodes[i] = ["C", [0,0,0]]
            origa.nodes[j] = ["C", [0,0,0]]
            origb.nodes[m] = ["C", [0,0,0]]
            origb.nodes[n] = ["C", [0,0,0]]
            if equalchem(origa, origb):
                x = nx.shortest_path_length(convert_to_networkx(a), source=i, target=j)
                y = nx.shortest_path_length(convert_to_networkx(b), source=m, target=n)
                if x == y:
                    continue
                if x < y:
                    return -1
                else:
                    return 1
                
                
    return 0
    
print("done")
print(len(done2))
done2 = sorted(done2, key=functools.cmp_to_key(cit))
print("done2")
for item in done2:
    image.append(display(item, False))
image = create_square_album(image)

image.show()
