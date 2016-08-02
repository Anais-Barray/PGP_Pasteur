#!/usr/bin/python
import sys,os

kaiju_classif_file=sys.argv[1]
kraken_labels_file=sys.argv[2]
search_space=open(kraken_labels_file,'r').read().splitlines()
output = open(kaiju_classif_file+".virus_classified", "w")
n=0

for line in open(kaiju_classif_file):
	if (line[0]=="C"):
		seq_id=line.split()[1]
		taxo=search_space[n].split("\t")[1]
		n+=1
		if ("Eukaryota" in taxo):
			newline="U"+line[1:]
			output.write(newline)
		elif ("Bacteria" in taxo):
			newline="U"+line[1:]
			output.write(newline)
		else:
			output.write(line)
	else:
		output.write(line)
output.close()
