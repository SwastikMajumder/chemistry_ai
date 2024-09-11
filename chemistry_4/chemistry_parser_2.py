import copy
from lark import Lark, Tree


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

# Define the grammar for simple arithmetic expressions
grammar = """
%import common.NUMBER
%import common.WS_INLINE
%ignore WS_INLINE

?start: expr

?expr: segment ("," segment)*

?segment: compound (";" compound)*

?compound: position ("," position)* "*" base
         | position ("," position)*
         
?base: VARIABLE | "(" segment ")"

?position: NUMBER

VARIABLE: "meth" | "eth" | "prop" | "but" | "pent" | "hex" | "cyclometh" | "cycloeth" | "cycloprop" | "cyclobut" | "cyclopent" | "cyclohex" | "hept" | "oct" | "non" | "dec" | "undec" | "dodec" | "tridec" | "tetradec" | "pentadec" | "hexadec" | "heptadec" | "octadec" | "nonadec" | "icos" | "cyclohept" | "cyclooct" | "cyclonon" | "cyclodec" | "cycloundec" | "cyclododec"

"""

# Create a parser
parser = Lark(grammar, start='start', parser='lalr')

# Example equation to parse
def take_input(equation):#"sin(3)+1+2"

  # Parse the equation
  parse_tree = parser.parse(equation)

  def convert_to_treenode(parse_tree):
      def tree_to_treenode(tree):
          if isinstance(tree, Tree):
              node = TreeNode(tree.data)
              node.children = [tree_to_treenode(child) for child in tree.children]
              return node
          else:  # Leaf node
              return TreeNode(str(tree))

      return tree_to_treenode(parse_tree)

  # Convert and print TreeNode structure
  tree_node = convert_to_treenode(parse_tree)
  tree_node = str_form(tree_node)
  return tree_node
  #print(tree_node)
  for item in ["add", "sin", "cos", "tan", "mul", "pow", "div"]:
    tree_node = tree_node.replace(str(item), "f_" + str(item))
  tree_node = tree_form(tree_node)
  for i in range(100,-1,-1):
    tree_node = replace(tree_node, tree_form(str(i)), tree_form("d_"+str(i)))
  for i in range(3):
    tree_node = replace(tree_node, tree_form(["x", "y", "z"][i]), tree_form("v_"+str(i)))
  tree_node = str_form(tree_node)
  return tree_node
