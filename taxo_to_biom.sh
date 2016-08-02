#!/bin/bash

### 	Ce programme transforme des fichiers de taxonomies (sorties Blast avec la taxonomie taxoptimizer) en format BIOM (JSON). 
###		Ce fichier peut ensuite etre analyse avec Phinch (version locale dans /Script/Phinch en lancant index.html), ou sur le site phinch.org.

#Prend en entree :
#1 - le fichier contenant la liste des fichiers (sortie taxoptimizer) a convertir en fichier BIOM [mandatory].
#  Il peut s'agir de fichiers de reads simples, ou de deux fichiers suite a un assemblage. 
#  Dans ce cas, mettre sur la meme ligne, separes par un ";", le fichier de contig suivi du fichier de reads non assemblés
	#exemple => list_file.txt: 
	#sample1:/path/to/file/file1.taxo
	#sample2:/path/to/file/file2.contigs.taxo;/path/to/file/file2.unassembled_reads.taxo
#2 - le fichier de metadata
	#exemple => map.tsv:
	#SampleID	Group   HealthCondition	SampleType	Date
	#sample1	1	severe	Blood 	2010-01-22
	#sample2	3	light	Blood   2010-02-16
	  #	/!\ le fichier de metadata doit suivre des normes, plus d'info sur le formattage de ce fichier aux liens suivants :
	  # 	https://github.com/PitchInteractiveInc/Phinch/wiki/Quick-Start
	  # 	https://github.com/PitchInteractiveInc/Phinch/wiki/Format
#3 - Le fichier de sortie, avec l'extension .biom

	
#Calling example : taxo_to_biom.sh -i comp_samples_otutable.tsv -o comp_samples_otutable.biom -m map.tsv


#	  .			THIS PROGRAM REQUIRES THAT BIOM HAS BEEN INSTALLED ON YOUR SESSION
#	 / \		Follow the instructions below as an example of installation on bic cluster:
#   / | \		 
#  /__°__\		module load Python/2.7.8	 	(or your favorite version of python with pip installed )
#				pip install --user numpy	 	(no need for --user option if installed on local)
#				pip install --user biom-format
#				pip install --user h5py
#
#
#	Program written by A.Barray
#	Date 27.07.2016

usage="Usage : $(basename "$0") [-h] [-i otutable_file] [-m metadata_file] [-o output_file] "

## PARSING OPTIONS
while getopts i:m:o: option
do
  case $option in
  i)	OTUTABLEFILE=$OPTARG
		if [ ! -e $OTUTABLEFILE ]; then echo "   problem with input file (option -i): '$OTUTABLEFILE' does not exist" ; exit ; fi;;

  m)	METADATAFILE=$OPTARG
		if [ ! -e $METADATAFILE ]; then echo "   problem with input file (option -m): '$METADATAFILE' does not exist" ; exit ; fi;;
		
  o)	OUTPUTFILE=$OPTARG
		oflag=true;;
		
  :)	printf "missing argument for -%s\n" "$OPTARG" >&2
		echo "$usage" >&2
		exit 1 ;;
  \?)	printf "illegal option: -%s\n" "$OPTARG" >&2
		echo "$usage" >&2
		exit 1 ;;
  h)	echo "$usage"
		exit 1;;
		
  esac
done
shift "$((OPTIND-1))"

if [[ -z $oflag ]]
then
	echo "Exit : Please indicate an output file" >&2
	echo "$usage";
	exit 1;
fi

# VARIABLES
OUTPUT_NAME="${OUTPUTFILE%.*}"
SCRIPTDIR=/pasteur/projets/specific/PGP_Work/Scripts

# TAXO TO OTU TABLE
echo "... Step 1 : Converting taxo files to OTU table ..."
qrsh -q cibu -cwd -now n -V ${SCRIPTDIR}/taxo_to_otutable.py -i ${OTUTABLEFILE} -o ${OUTPUT_NAME}.tsv
echo
# OTU TABLE TO BIOM
echo "... Step 2 : Converting OTU table to BIOM file ..."
qrsh -q cibu -cwd -now n -V ${SCRIPTDIR}/otutable_to_biom.sh -i ${OUTPUT_NAME}.tsv -m ${METADATAFILE} -o ${OUTPUT_NAME}.biom

