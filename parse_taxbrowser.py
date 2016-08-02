 #!/usr/bin/python

import os
import sys

# Definition of the classe Node
class Node:
	"""Noeud"""
	def __init__(self):
		self.tax_id = 0       # Number of the tax id.
		self.parent = 0       # Number of the parent of this node
		self.children = []    # List of the children of this node
		self.tip = 0          # Tip=1 if it's a terminal node, 0 if not.
		self.name = ""        # Name of the node: taxa if it's a terminal node, numero if not.      
	def genealogy(self):      # Trace genealogy from root to leaf
		ancestors = []        # Initialise the list of all nodes from root to leaf.
		tax_id = self.tax_id  # Define leaf
		while 1:
			if name_object.has_key(tax_id):
				ancestors.append(tax_id)
				tax_id = name_object[tax_id].parent
			else:
				break
			if tax_id == "1":
				# If it is the root, we reached the end.
				# Add it to the list and break the loop
				ancestors.append(tax_id)
				break
		return ancestors # Return the list

# Function to find common ancestor between two nodes or more 

def common_ancestor(node_list):
	global name_object
	list1 = name_object[node_list[0]].genealogy()  # Define the whole genealogy of the first node
	for node in node_list:
		list2 = name_object[node].genealogy()      # Define the whole genealogy of the second node
		ancestral_list = []                            
		for i in list1:
			if i in list2:                         # Identify common nodes between the two genealogy
				ancestral_list.append(i)                
		list1 = ancestral_list                     # Reassing ancestral_list to list 1.
	common_ancestor = ancestral_list[0]            # Finally, the first node of the ancestra_list is the common ancestor of all nodes.
	return common_ancestor                         # Return a node 

def all_descendant(self): # Find all children from a node
		terminal_nodes  = []
		list_descendant = []
		list_descendant.append(self.tax_id)
		for i in list_descendant:
			#print i
			if name_object[i].tip == 1:
				terminal_nodes.append(i)
			else:
				for j in name_object[i].children:
					list_descendant.append(j)
		return terminal_nodes # Return a list 



#############################
#                           #
#   Read taxonomy files     #
#                           #
#############################

######################
#
# Load names defintion

name_dict = {}          # Initialise dictionary with TAX_ID:NAME
#~ name_dict_reverse = {}  # Initialise dictionary with NAME:TAX_ID

# Load  NCBI names file ("names.dmp")
name_file =  open("names.dmp","r")
while 1:
	line = name_file.readline()
	if line == "":
		break
	line = line.rstrip()
	line = line.replace("\t","")
	tab = line.split("|")
	if tab[3] == "scientific name":
		tax_id, name = tab[0], tab[1]     # Assign tax_id and name ...
		name_dict[tax_id] = name          # ... and load them
		#name_dict_reverse[name] = tax_id  # ... into dictionaries
name_file.close()


######################
#
# Load taxonomy

# Define taxonomy variable
global name_object
name_object = {}


# Load taxonomy NCBI file ("nodes.dmp")
taxonomy_file = open("nodes.dmp","r")
while 1:
	line = taxonomy_file.readline()
	if line == "":
		break
	#print line
	line = line.replace("\t","")
	tab = line.split("|")
   
	tax_id = str(tab[0])
	tax_id_parent = str(tab[1])
	division = str(tab[4])

	# Define name of the taxid
	name = "unknown"
	if tax_id in name_dict:
		name = name_dict[tax_id]
   
	if not name_object.has_key(tax_id):
		name_object[tax_id] = Node()
	name_object[tax_id].tax_id   = tax_id        # Assign tax_id
	name_object[tax_id].parent   = tax_id_parent # Assign tax_id parent
	name_object[tax_id].name     = name          # Assign name
   
	if  tax_id_parent in name_object:
		children = name_object[tax_id_parent].children  # If parent is is already in the object
		children.append(tax_id)                  # ...we found its children.
		name_object[tax_id_parent].children = children  # ... so add them to the parent
taxonomy_file.close()

##### Print taxo #####
for name_id in name_dict :
	print name_id, " : ", name_dict[name_id]
