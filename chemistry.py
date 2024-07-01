class Graph:
    def __init__(self):
        self.nodes = []
        self.edges = []
    
    def add_node(self, node_id):
        if node_id not in self.nodes:
            self.nodes.append(node_id)
    
    def add_edge(self, node1_id, node2_id, value=None):
        self.edges.append([[node1_id, node2_id], value])

def terminal_carbon(compound):
    output = []
    for carbon in compound.nodes:
        adjacent_carbon = 0
        for bond in compound.edges:
            if carbon in bond[0]:
                adjacent_carbon += 1
        if adjacent_carbon == 1:
            output.append(carbon)
    return output

all_chain_output = []
sub_chain_output = []
def all_chain(chain, compound, sub_chain):
    end_reached = True
    atom_iter = []
    for atom in compound.nodes:
        if atom not in chain and {atom, chain[-1]} in [set(item[0]) for item in compound.edges]:
            atom_iter.append(atom)
    if len(atom_iter) == 0:
        all_chain_output.append(chain)
        sub_chain_output.append(sub_chain)
    elif len(atom_iter) == 1:
        atom = atom_iter[0]
        all_chain(chain + [atom], compound, sub_chain)
    else:
        for atom in atom_iter:
            all_chain(chain + [atom], compound, sub_chain + [len(chain)]*(len(atom_iter)-1))
        

graph = Graph()

# Add nodes
graph.add_node(1)
graph.add_node(2)
graph.add_node(3)
graph.add_node(4)
graph.add_node(5)
graph.add_node(6)

# Add edges with values
graph.add_edge(1, 2, value=1)
graph.add_edge(2, 3, value=1)
graph.add_edge(3, 4, value=1)
graph.add_edge(3, 5, value=1)
graph.add_edge(5, 6, value=1)

for atom in terminal_carbon(graph):
    all_chain([atom], graph, [])
for i in range(len(all_chain_output)):
    print(all_chain_output[i], sub_chain_output[i])
