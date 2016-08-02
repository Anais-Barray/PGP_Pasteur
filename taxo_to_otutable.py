#!/usr/bin/env python3

#Ce programme transforme des fichiers de taxonomies (sorties Blast avec la taxonomie taxoptimizer) en OTU table (TSV). 
#Etape necessaire necessaire pour produire un fichier au format BIOM avec le programme 'otutable_to_biom.sh', qui peut ensuite etre analyse avec Phinch.

#Prend en entree :
#1 - le fichier contenant la liste des fichiers (sortie taxoptimizer) a convertir en fichier BIOM [mandatory].
#  Il peut s'agir de fichiers de reads simples, ou de deux fichiers suite a un assemblage. 
#  Dans ce cas, mettre sur la meme ligne, separes par un ";", le fichier de contig suivi du fichier de reads non assemblés
	#exemple => list_file.txt: 
	#sample1:/path/to/file/file1.taxo
	#sample2:/path/to/file/file2.contigs.taxo;/path/to/file/file2.unassembled_reads.taxo
#2 - le fichier de sortie [optional].
	#exemple => comparison_samples_OTUtable.tsv

#Calling example : python3 taxo_to_otutable.py -i list_sample.txt -o comparison_samples_OTUtable.tsv

#	Program written by A.Barray
#	Date 25.07.2016

import os
import sys
import string
import argparse
import re


################################################################################################
####################################	FUNCTIONS	 ###########################################
################################################################################################

def getSampleNames(F):
	# returns a list of samples' names, written before ":" in input.
	
	table= []
	fofs=open(F,'r')
	
	while(True):
		line=fofs.readline()
		if not line: break	#EOF
		if not line.strip(): continue # avoid empty lines
		table.append(line.split(":")[0].strip())
	fofs.close()
	
	return table
	
	
def getDicoFiles(F):
	# returns a dictionnary of samples, where the keys are sample names, and values are path/to/file.

	dicoFiles = {}
	fofs=open(F,'r')
	
	while(True):
		line=fofs.readline()
		if not line: break	#EOF
		if not line.strip(): continue # avoid empty lines
		sample_name=line.split(":")[0]
		paths=line.split(":")[1]
		tab_line=paths[:-1].split(";") # [:-1] removes the \n
		dicoFiles[sample_name]=tab_line
	fofs.close()
	
	return dicoFiles


def getDicoOTUtotaxo(dicoFiles):	
	# returns a dictionnary containing the formatted taxonomies encountered in all samples. 
	# to each taxonomy (value) is attributed an arbitrary OTU number (key).
	# ex : { 1 : "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__Salmonella; s__Salmonella enterica" ; 
	#		 2 : "k__Viruses; p__; c__; o__; f__Bunyaviridae; g__Orthobunyavirus; s__Tete virus" ... } 
	
	dicoOTUtotaxo={}
	dico_species={}
	dico_species_taxo={}
	OTUcount=1
	
	for sample_name in sorted(dicoFiles.keys()):
		for sample_file_path in dicoFiles[sample_name] :
			with open(sample_file_path,'r') as f :
				
				lines=f.readlines()
				for line in lines:
					line=line.rstrip().split("\t")
					if (len(line)<13) : # no taxonomy fields
						continue
					else:
						taxo_line=line[13]
						species=line[12]
						
						if "(species)" in taxo_line :
							s=re.search(".*; (.*) \(species\)",taxo_line)
							sp=s.group(1)
						else :
							sp=species
							
						try :
							dico_species[sp] == 1 # species already listed in the dictionary of species
						except KeyError :
							dico_species[sp] = 1 # new species !
		
							formated_taxo=""	# formatted taxo for .biom files : k__; p__; c__; o__; f__; g__; s__
							
							if "Viruses" in taxo_line :
								formated_taxo=formated_taxo+"k__Viruses; "
							elif "(superkingdom)" in taxo_line :
								k=re.search("; (.*) \(superkingdom\)",taxo_line)
								formated_taxo="k__"+k.group(1)+"; "
							else:
								formated_taxo="k__; "
							
							if "(phylum)" in taxo_line :
								p=re.search(".*; (.*) \(phylum\)",taxo_line)
								formated_taxo=formated_taxo+"p__"+p.group(1)+"; "
							else:
								formated_taxo=formated_taxo+"p__; "
		
							if "(class)" in taxo_line :
								c=re.search(".*; (.*) \(class\)",taxo_line)
								formated_taxo=formated_taxo+"c__"+c.group(1)+"; "
							else:
								formated_taxo=formated_taxo+"c__; "
								
							if "(order)" in taxo_line :
								o=re.search(".*; (.*) \(order\)",taxo_line)
								formated_taxo=formated_taxo+"o__"+o.group(1)+"; "
							else:
								formated_taxo=formated_taxo+"o__; "
								
							if "(family)" in taxo_line :
								f=re.search(".*; (.*) \(family\)",taxo_line)
								formated_taxo=formated_taxo+"f__"+f.group(1)+"; "
							else:
								formated_taxo=formated_taxo+"f__; "
		
							if "(genus)" in taxo_line :
								g=re.search(".*; (.*) \(genus\)",taxo_line)
								formated_taxo=formated_taxo+"g__"+g.group(1)+"; "
							else:
								formated_taxo=formated_taxo+"g__; "
		
							formated_taxo=formated_taxo+"s__"+sp
									
							dicoOTUtotaxo[OTUcount] = [sp,formated_taxo]
							OTUcount += 1
						species=""
						sp=""

	return dicoOTUtotaxo
	

def getDicoOTUcount(sample_list_files, dicoOTUtotaxo):
	# returns a dictionnary containing for each taxonomy (key) the occurence (value, has to be a float), i.e the number of times this taxonomy is encountered in a sample
	
	dicoOTUcount={}
	count=0
	
	if (len(sample_list_files)>1): 
		contig_file=sample_list_files[0]
		unassembled_file=sample_list_files[1]
	
		for OTU in sorted(dicoOTUtotaxo.keys()):
			sp=dicoOTUtotaxo[OTU][0]
			with open(contig_file,'r') as f :
				lines=f.readlines()
				for line in lines:
					line_tab=line.rstrip().split("\t")
					if sp==line_tab[12]:
						nb=re.search(".*_nbreads_(.*)_cov.*", line)
						count+=int(nb.group(1))
					elif sp in line_tab[13]:
						nb=re.search(".*_nbreads_(.*)_cov.*", line)
						count+=int(nb.group(1))
			with open(unassembled_file,'r') as g :
				lines=g.readlines()
				for line in lines:
					line_tab=line.rstrip().split("\t")
					if sp==line_tab[12]:
						count+=1
					elif sp in line_tab[13]:
						count+=1
			dicoOTUcount[OTU]="%.1f" % count
			count=0
		
	else : 
		sample_file=sample_list_files[0]

		for OTU in sorted(dicoOTUtotaxo.keys()):
			sp=dicoOTUtotaxo[OTU][0]
			with open(sample_file,'r') as f :
				lines=f.readlines()
				for line in lines:
					line_tab=line.rstrip().split("\t")
					if sp==line_tab[12]:
						count+=1
					elif sp in line_tab[13]:
						count+=1
			dicoOTUcount[OTU]="%.1f" % count
			count=0
		
	return dicoOTUcount

################################################################################################
########################################	MAIN	 ###########################################
################################################################################################

def main():
	
	# Parsing options
	parser = argparse.ArgumentParser(description='Converts a list of taxonomy files (Blast + Taxoptimizer output) to a .biom format', add_help=True)
	
	parser.add_argument("-i", type=str, dest="infile", required=True,
						help="input file of files (a line format type = sample_name:/path/to/file/filename)" )

	parser.add_argument("-o", type=str, dest="outfile",
						help="output file as an OTU table (TSV) format.", default="output_taxo_to_otu_table.tsv")
	
	#TOADD : 
	#	-min value / max value or % to retain taxonomies. 
	
	args = parser.parse_args()
	
	# Variables
	input_file=str(args.infile)
	output_file=str(args.outfile)
	
	# Generating tables and dictionaries
	print("... Producing preliminary tables and dictionaries ...")
	list_sample_names=getSampleNames(input_file)
	dico_sample=getDicoFiles(input_file)
	dico_OTU_to_taxo=getDicoOTUtotaxo(dico_sample)
	

	print("... Creating matrix of OTU counts for all samples ...")
	dico_allsample_OTUcount={}
	for sample in list_sample_names :
		print("Filling with", sample)
		dico_OTU_count=getDicoOTUcount(dico_sample[sample],dico_OTU_to_taxo)
		dico_allsample_OTUcount[sample]=dico_OTU_count
	
	#Generating output
	print("... Generating Output ...")
	ofile=open(output_file, "w")
	first_line="#OTU ID\t"+"\t".join(sorted(list_sample_names))+"\ttaxonomy\n"
	ofile.write(first_line)
	for OTU in sorted(dico_OTU_to_taxo.keys()):
		otu_line=str(OTU)
		sum_otu=0
		for sample in sorted(dico_allsample_OTUcount.keys()):
			count = float(dico_allsample_OTUcount[sample][OTU])
			sum_otu+=count
			otu_line=otu_line+"\t"+str(count)
		otu_line=otu_line+"\t"+dico_OTU_to_taxo[OTU][1]+"\n"
		if (sum_otu >= 10) :
			ofile.write(otu_line)

	ofile.close()
	
	
if __name__ == "__main__":
	main()
