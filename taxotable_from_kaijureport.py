#!/usr/bin/python3
import os,sys,re

############################ FUNCTIONS #################################

def taxo_dict_per_sample(sample):
	taxo_dict={}
	sample1=sample
	
	with open(sample1,'r') as f :
		lines=f.readlines()
		lines=lines[2:-3]
		for line in lines:
			taxo=line.rstrip().split("\t")[2]
			nb_read=line.rstrip().split()[1]
			taxo_dict[taxo]=int(nb_read)
	
	return taxo_dict

########################################################################
############################## MAIN #################################### 
########################################################################

def main():
	directory = '/pasteur/projets/specific/PGP_Work/Analyses/138_ANSM_VSauvage'
	extension = 'virus.genus-report'
	sample_list=[]
	sample_dict={}
	all_list_taxo=[]
	all_sample_dict={}
	output=open("taxo_table_from_kaiju_report_virus.txt",'w')
	
	#Create a dictionary of samples
	for dirpath, dirnames, files in os.walk(directory):
		for name in files:
			if name.lower().endswith(extension):
				sample_list.append(os.path.join(dirpath, name))
				
	   
	#List of all taxa
	for sample_name in sample_list:
		taxo_dict=taxo_dict_per_sample(sample_name)
		list_taxo=taxo_dict.keys()
		for taxo in list_taxo:
			if not taxo in all_list_taxo:
				all_list_taxo.append(taxo)
		
	#Fill the dictionary containing, for each sample, each abundance of each taxonomy.
	for sample_name in sample_list:
		all_sample_dict[sample_name]={}
		taxo_dict=taxo_dict_per_sample(sample_name)
		
		for taxo in sorted(all_list_taxo):
			if taxo in taxo_dict.keys():
				all_sample_dict[sample_name][taxo]=taxo_dict.get(taxo)
			else:
				all_sample_dict[sample_name][taxo]=0
				
	#### Writing in output file
	#first line = sample names
	taxo_line="\t".join(sorted(all_list_taxo))
	taxo_line="\t"+taxo_line+"\n"
	output.write(taxo_line)
	
	#next lines = taxo abundance per sample
	for sample in sorted(all_sample_dict):
		sample_line=sample+"\t"
		for taxo in sorted(all_list_taxo):
			sample_line+=str(all_sample_dict[sample][taxo])+"\t"
		sample_line+="\n"
		output.write(sample_line)
	
	

		#~ for sample in sorted(all_sample_dict):
			#~ for taxo in sorted(all_list_taxo):
				#~ taxo_line=taxo+"\t"
			#~ taxo_line+=str(all_sample_dict[sample][taxo])+"\t"
		#~ taxo_line+="\n"
		#~ output.write(taxo_line)
	#~ 
	
	output.close()
	
	sys.exit(1)

########################################################################

if __name__ == "__main__":
	main()
		
