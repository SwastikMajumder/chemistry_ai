import copy

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

def is_number(s):
    try:
        float(s)  # Try to convert the string to a float
        return True
    except ValueError:
        return False

# Generate transformations of a given equation provided only one formula to do so
# We can call this function multiple times with different formulas, in case we want to use more than one
# This function is also responsible for computing arithmetic, pass do_only_arithmetic as True (others param it would ignore), to do so
def apply_individual_formula_on_given_equation(equation, formula_lhs, formula_rhs, do_only_arithmetic=False, structure_satisfy=False):
    variable_list = {}
    
    def node_type(s):
        if s[:2] == "f_":
            return s
        else:
            return s[:2]
    def does_given_equation_satisfy_forumla_lhs_structure(equation, formula_lhs):
        nonlocal variable_list
        # u can accept anything and p is expecting only integers
        # if there is variable in the formula
        if node_type(formula_lhs.name) in {"u_", "p_"}: 
            if formula_lhs.name in variable_list.keys(): # check if that variable has previously appeared or not
                return str_form(variable_list[formula_lhs.name]) == str_form(equation) # if yes, then the contents should be same
            else: # otherwise, extract the data from the given equation
                if node_type(formula_lhs.name) == "p_" and "v_" in str_form(equation): # if formula has a p type variable, it only accepts integers
                    return False
                variable_list[formula_lhs.name] = copy.deepcopy(equation)
                return True
        if equation.name != formula_lhs.name or len(equation.children) != len(formula_lhs.children): # the formula structure should match with given equation
            return False
        for i in range(len(equation.children)): # go through every children and explore the whole formula / equation
            if does_given_equation_satisfy_forumla_lhs_structure(equation.children[i], formula_lhs.children[i]) is False:
                return False
        return True
    if structure_satisfy:
      return does_given_equation_satisfy_forumla_lhs_structure(equation, formula_lhs)
    # transform the equation as a whole aka perform the transformation operation on the entire thing and not only on a certain part of the equation
    def formula_apply_root(formula):
        nonlocal variable_list
        if formula.name in variable_list.keys():
            return variable_list[formula.name] # fill the extracted data on the formula rhs structure
        data_to_return = TreeNode(formula.name, None) # produce nodes for the new transformed equation
        for child in formula.children:
            data_to_return.children.append(formula_apply_root(copy.deepcopy(child))) # slowly build the transformed equation
        return data_to_return
    count_target_node = 1
    # try applying formula on various parts of the equation
    def formula_apply_various_sub_equation(equation, formula_lhs, formula_rhs, do_only_arithmetic):
        nonlocal variable_list
        nonlocal count_target_node
        data_to_return = TreeNode(equation.name, children=[])
        variable_list = {}
        if do_only_arithmetic == False:
            if does_given_equation_satisfy_forumla_lhs_structure(equation, copy.deepcopy(formula_lhs)) is True: # if formula lhs structure is satisfied by the equation given
                count_target_node -= 1
                if count_target_node == 0: # and its the location we want to do the transformation on
                    return formula_apply_root(copy.deepcopy(formula_rhs)) # transform
        else: # perform arithmetic
            if len(equation.children) == 2 and all(node_type(item.name) == "d_" and item.name[2:]for item in equation.children): # if only numbers
                x = []
                for item in equation.children:
                    x.append(float(item.name[2:])) # convert string into a number
                if equation.name == "f_add":
                    count_target_node -= 1
                    if count_target_node == 0: # if its the location we want to perform arithmetic on
                        return TreeNode("d_" + str(sum(x))) # add all
                elif equation.name == "f_mul":
                    count_target_node -= 1
                    if count_target_node == 0:
                        p = 1
                        for item in x:
                            p *= item # multiply all
                        return TreeNode("d_" + str(p))
                elif equation.name == "f_pow": # power should be two or a natural number more than two
                    count_target_node -= 1
                    if count_target_node == 0:
                        return TreeNode("d_"+str(float(x[0]**x[1])))
                elif equation.name == "f_sub":
                    count_target_node -= 1
                    if count_target_node == 0: # if its the location we want to perform arithmetic on
                        return TreeNode("d_" + str(x[1]-x[0]))
                elif equation.name == "f_div":
                    count_target_node -= 1
                    if count_target_node == 0: # if its the location we want to perform arithmetic on
                        return TreeNode("d_" + str(float(x[0]/x[1])))
        if node_type(equation.name) in {"d_", "v_"}: # reached a leaf node
            return equation
        for child in equation.children: # slowly build the transformed equation
            data_to_return.children.append(formula_apply_various_sub_equation(copy.deepcopy(child), formula_lhs, formula_rhs, do_only_arithmetic))
        return data_to_return
    cn = 0
    # count how many locations are present in the given equation
    def count_nodes(equation):
        nonlocal cn
        cn += 1
        for child in equation.children:
            count_nodes(child)
    transformed_equation_list = []
    count_nodes(equation)
    for i in range(1, cn + 1): # iterate over all location in the equation tree
        count_target_node = i
        orig_len = len(transformed_equation_list)
        tmp = formula_apply_various_sub_equation(equation, formula_lhs, formula_rhs, do_only_arithmetic)
        if str_form(tmp) != str_form(equation): # don't produce duplication, or don't if nothing changed because of transformation impossbility in that location
            transformed_equation_list.append(str_form(tmp)) # add this transformation to our list
    return transformed_equation_list 

# Function to generate neighbor equations
def generate_transformation(equation, file_name):
    input_f, output_f = return_formula_file(file_name) # load formula file
    transformed_equation_list = []
    for i in range(len(input_f)): # go through all formulas and collect if they can possibly transform
        transformed_equation_list += apply_individual_formula_on_given_equation(tree_form(copy.deepcopy(equation)), copy.deepcopy(input_f[i]), copy.deepcopy(output_f[i]))
    return list(set(transformed_equation_list)) # set list to remove duplications

# Function to generate neighbor equations
def generate_arithmetical_transformation(equation):
    transformed_equation_list = []
    transformed_equation_list += apply_individual_formula_on_given_equation(tree_form(equation), None, None, True) # perform arithmetic
    return list(set(transformed_equation_list)) # set list to remove duplications

# Function to read formula file
def return_formula_file(file_name):
    global add_eq
    with open(file_name, 'r') as file:
      content = file.read()
    x = content.split("\n\n")
    #x += add_eq.split("\n\n")
    #print("\n\n".join(x))
    input_f = [x[i] for i in range(0, len(x), 2)] # alternative formula lhs and then formula rhs
    output_f = [x[i] for i in range(1, len(x), 2)]
    input_f = [tree_form(item) for item in input_f] # convert into tree form
    output_f = [tree_form(item) for item in output_f]
    return [input_f, output_f] # return

def search(equation, depth, file_list, auto_arithmetic=True, visited=None):
    if depth == 0: # limit the search
        return None
    if visited is None:
        visited = set()

    print(string_equation(equation))
    #print(equation)
    if equation in visited:
        return None
    visited.add(equation)
    output =[]
    if file_list[0]:
      output += generate_transformation(equation, file_list[0])
    if auto_arithmetic:
      output += generate_arithmetical_transformation(equation)
    if len(output) > 0:
      output = [output[0]]
    else:
      if file_list[1]:
        output += generate_transformation(equation, file_list[1])
      if not auto_arithmetic:
        output += generate_arithmetical_transformation(equation)
      if file_list[2] and len(output) == 0:
          output += generate_transformation(equation, file_list[2])
    for i in range(len(output)):
        result = search(output[i], depth-1, file_list, auto_arithmetic, visited) # recursively find even more equals
        if result is not None:
            output += result # hoard them
    output = list(set(output))
    return output

sign = {
    "f_add": "+", 
    "f_mul": "*", 
    "f_pow": "^", 
    "f_div": "/", 
    "f_eq": "=", 
    "f_amu": "amu"
}

# Helper function to convert the equation tree to a string ensuring binary operations
def string_equation_helper(equation_tree):
    if equation_tree.children == []:
        return equation_tree.name  # leaf node
    
    operation = equation_tree.name  # get the operation name
    operator_symbol = sign.get(operation, operation[2:])  # get the symbol or default to the operation name

    # Assuming it's binary, take the first two children
    if len(equation_tree.children) == 2:
        left_operand = string_equation_helper(copy.deepcopy(equation_tree.children[0]))
        right_operand = string_equation_helper(copy.deepcopy(equation_tree.children[1]))
        return f"({left_operand} {operator_symbol} {right_operand})"
    
    # If there are more than two children, combine them pairwise
    s = string_equation_helper(copy.deepcopy(equation_tree.children[0]))
    for i in range(1, len(equation_tree.children)):
        s = f"({s} {operator_symbol} {string_equation_helper(copy.deepcopy(equation_tree.children[i]))})"
    
    return s

# Main function to convert the equation to its string representation
def string_equation(eq): 
    eq = eq.replace("v_0", "x")  # replace v_ with #
    eq = eq.replace("d_", "")   # remove d_ prefix
    eq = eq.replace("c_", "")
    
    return string_equation_helper(tree_form(eq))

def replace(equation, find, r):
  if str_form(equation) == str_form(find):
    return r
  col = TreeNode(equation.name, [])
  for child in equation.children:
    col.children.append(replace(child, find, r))
  return col

import numpy as np

def convert_reaction(reaction, compound_data):
    # Split the reaction into reactants and products
    reaction = reaction.replace(" ", "")
    reactants, products = reaction.split("->")
    
    # Convert reactants and products into lists of dictionaries
    reactants_list = [compound_data[compound] for compound in reactants.split("+")]
    products_list = [compound_data[compound] for compound in products.split("+")]
    
    # Combine reactants and products into the desired format
    return [reactants_list, products_list]

def convert_to_string(reaction_data, x, compound_data):
    x = x[0]+x[1]
    # Create a reverse mapping from element composition to compound name
    reverse_compound_data = {frozenset(v.items()): k for k, v in compound_data.items()}
    
    # Extract reactants and products with their coefficients
    reactants, products = reaction_data
    reactant_coeffs = x[:len(reactants)]
    product_coeffs = x[len(reactants):]
    
    # Build the string for reactants
    reactant_strs = [f"{coeff} {reverse_compound_data[frozenset(compound.items())]}"
                     for coeff, compound in zip(reactant_coeffs, reactants)]
    
    # Build the string for products
    product_strs = [f"{coeff} {reverse_compound_data[frozenset(compound.items())]}"
                    for coeff, compound in zip(product_coeffs, products)]
    
    # Combine into the final reaction string
    reaction_string = " + ".join(reactant_strs) + " -> " + " + ".join(product_strs)
    
    return reaction_string.replace("1 ", "")

# Example usage
compound_data = {
    "C3H8": {"C": 3, "H": 8},
    "O2": {"O": 2},
    "CO2": {"C": 1, "O": 2},
    "H2O": {"H": 2, "O": 1},
    "MgO": {"Mg": 1, "O": 1},
    "Mg": {"Mg": 1}
}

#reaction = [[{"C": 3, "H": 8}, {"O": 2}], [{"C": 1, "O": 2}, {"H": 2, "O": 1}]]
reaction = "C3H8+O2->CO2+H2O"

periodic_table = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg"]
atomic_mass = [1, 4, 7, 9, 10.8, 12, 14, 18, 19, 20.2, 23, 24.3]

def molar_mass(compound):
    mass = 0
    for item in compound.keys():
        mass += atomic_mass[periodic_table.index(item)]*compound[item]
    return mass

def convert2gram(quantity):
    if quantity[-1] == "mol":
        return [quantity[0], quantity[1]*molar_mass(quantity[0]), quantity[-1]]
    elif quantity[-1] == "gram":
        return quantity
    elif quantity[-1] == "kilogram":
        return [quantity[0], quantity[1]*1000, quantity[-1]]
    return quantity

def convert2mol(quantity):
    if quantity[-1] == "gram":
        return [quantity[0], quantity[1]/molar_mass(quantity[0]), quantity[-1]]
    elif quantity[-1] == "mol":
        return quantity
    elif quantity[-1] == "kilogram":
        return [quantity[0], quantity[1]/1000, quantity[-1]]
    return quantity

def balance(reaction):
  
  element_present = [[],[]]
  for i in range(2):
    for item in reaction[i]:
      element_present += list(item.keys())
  matrix = []
  for element in periodic_table:
    if element not in element_present:
      continue
    matrix.append([0]*(len(reaction[0])+len(reaction[1])))
    for i in range(2):
      for index, compound in enumerate(reaction[i]):
        if element not in compound.keys():
          continue  
        if i == 0:
          matrix[-1][index] = compound[element]
        else:
          matrix[-1][index+len(reaction[0])] = -compound[element]
  
  B = [0]*len(matrix)+[1]
  matrix.append([0]*(len(reaction[0])+len(reaction[1])))
  matrix[-1][0] = 1
  A = np.array(matrix)
  B = np.array(B)
  output = [x for x in np.linalg.solve(A, B)]
  output = [output[:len(reaction[0])], output[len(reaction[0]):]]
  #print(output)
  return output

def find_quantity(reaction, coefficient, quantity_list):
    orig = quantity_list
    for j in range(len(quantity_list[1])):
        if None == quantity_list[1][j]:
            val = 999
            best = None
            for i in range(len(quantity_list[0])):
                if None != quantity_list[0][i]:
                    quantity = convert2mol(quantity_list[0][i])
                    quantity = [reaction[1][j], (quantity[1]/coefficient[0][i])*coefficient[1][j], "mol"]
                    if quantity[1] < val:
                        best = quantity
                        val = quantity[1]
            quantity_list[1][j] = best
    #if all(x==None for x in quantity_list[0]):
    for i in range(len(quantity_list[0])):
        if None == quantity_list[0][i]:
            for j in range(len(quantity_list[1])):
                if quantity_list[1][j] != None:
                    quantity = convert2mol(quantity_list[1][j])
                    quantity_list[0][i] = [reaction[0][i], (quantity[1]/coefficient[1][j])*coefficient[0][i], "mol"]
    #print(quantity_list)
    return quantity_list


reaction = convert_reaction("Mg + O2 -> MgO",compound_data)
q = [[[compound_data["Mg"], 2.4, "gram"], [compound_data["O2"], 3.2, "gram"]], [None]]
#q = find_quantity([[{"C": 3, "H": 8}, {"O": 2}], [{"C": 1, "O": 2}, {"H": 2, "O": 1}]], [[1,5],[3,4]], q)
eq = balance(reaction)
q = find_quantity(reaction, eq, q)
q = [[convert2gram(x) for x in q[0]],[convert2gram(x) for x in q[1]]]
print("Mg + O2 -> MgO")
print(convert_to_string(reaction, eq, compound_data))
print(q)

def flatten_tree(node):
    if not node.children:
        return node
    if node.name in {"f_add", "f_mul"}: # commutative property supporting functions
        merged_children = [] # merge all the children
        for child in node.children:
            flattened_child = flatten_tree(child)
            if flattened_child.name == node.name:
                merged_children.extend(flattened_child.children)
            else:
                merged_children.append(flattened_child)
        return TreeNode(node.name, merged_children)
    else:
        node.children = [flatten_tree(child) for child in node.children]
        return node
    
def reduce_single_equation(eq_list):
    output = []
    for eq in eq_list:
        print(string_equation(eq))
        eq2 = search(eq, 4, [None, "formula-list/single_equation.txt", None])
        for item in eq2 + [eq]:
            tmp = is_linear_equation(tree_form(item))
            if tmp is not None:
                output.append(tmp)
    all_key = []
    for item in output:
        #print(output)
        all_key += item.keys()
    all_key = list(set(all_key) - set(["const"]))
    matrix = []
    answer = []
    for i in range(len(output)):
        matrix.append([0]*len(all_key))
    for i in range(len(all_key)):
        for j in range(len(output)):
            if all_key[i] in output[j].keys():
                matrix[j][i] = float(output[j][all_key[i]][2:])
    for j in range(len(output)):
        answer.append(-float(output[j]["const"][2:]))
    result, residuals, rank, s = np.linalg.lstsq(np.array(matrix), np.array(answer), rcond=None)
    for i in range(len(all_key)):
        print((all_key[i] + " = " + str(list(result)[i])).replace("v_", ""))
    #print(result)
    #print(all_key)
    #
def solve_1(eq_list):
    for item in eq_list:
        tmp = reduce_single_equation(item)
def is_linear_equation(eq):
    if eq.name != "f_eq":
        return None
    if eq.children[1].name[:2] != "d_":
        return None
    eq = flatten_tree(eq)
    if eq.children[0].name[:2] == "v_":
        return {eq.children[0].name: "d_1", "const": eq.children[1].name}
    if eq.children[0].name == "f_mul" and eq.children[0].children[0].name[:2] == "d_" and eq.children[0].children[1].name[:2] == "v_":
        return {eq.children[0].children[1].name: eq.children[0].children[0].name, "const": eq.children[1].name}
    if eq.children[0].name != "f_add":
        return None
    def is_term(eq):
        if eq.name[:2] == "d_":
            return "const", eq.name
        if eq.name[:2] == "v_":
            return eq.name, "d_1"
        if eq.name != "f_mul":
            return None
        if len(eq.children) != 2:
            return None
        if eq.children[0].name[:2] != "d_":
            return None
        if eq.children[1].name[:2] != "v_":
            return None
        else:
            if str_form(eq).count(eq.children[1].name) != 1:
                return None
        return eq.children[1].name, eq.children[0].name
    output = {"const": eq.children[1].name}
    for child in eq.children[0].children:
        tmp2 = is_term(child)
        if tmp2 is not None:
            output[tmp2[0]] = tmp2[1]
        else:
            return None
    return output

eq_2 = """f_eq
 f_mul
  d_1.2
  v_0
 d_100"""
#reduce_single_equation(eq_2)
#tmp = convert_reaction(reaction,compound_data)
#eq = balance(tmp)
#print(tmp)
#print(convert_to_string(tmp, eq, compound_data))

eq = """f_mix
 v_0
 v_1
 v_2
 v_3
 v_4
 v_5
 v_6
 0.2
 d_1.25
 f_compound
  v_7
  v_8
  v_9
 f_compound
  v_10
  v_11
  v_12"""
print("a 20% m/v solution has density 1.25g/ml. tell %m/m for the solution")
print(string_equation(eq))
eq = """f_mix
 d_100
 volume_of_solution
 d_0.2
 molality_of_solution
 mole_fraction_of_solute
 mole_fraction_of_solvent
 weight_by_weight_ratio
 weight_by_volume_ratio
 d_1.2
 f_compound
  d_KCl
  mass_of_solute
  mole_of_solute
 f_compound
  d_H2O
  mass_of_solvent
  mole_of_solvent"""
print()
print("find the number of moles of KCl in a 100g 0.2 molar KCl solution")
print()
eq = """f_mix
 d_100
 v_volume_of_solution
 d_0.2
 v_molality_of_solution
 v_mole_fraction_of_solute
 v_mole_fraction_of_solvent
 v_weight_by_weight_ratio
 v_weight_by_volume_ratio
 d_1200
 f_compound
  c_KCl
  v_mass_of_solute
  v_mole_of_solute
 f_compound
  c_H2O
  v_mass_of_solvent
  v_mole_of_solvent"""
#print(string_equation(eq))
eq = sorted(search(eq, 10, [None, "formula-list/medium-main.txt", None]), key=lambda x: len(x))
reduce_single_equation(eq)
#print(string_equation(eq))

eq = """f_eq
 d_0.2
 f_div
  v_0
  d_0.1"""
search(eq, 4, [None, "formula-list/single_equation_3.txt", None])


