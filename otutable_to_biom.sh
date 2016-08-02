#!/bin/bash

### 	Programme convertissant une table d'OTU en format BIOM (JSON)
###		Ce fichier peut ensuite etre analyse avec Phinch

#Prend en entree :
#1 - Le fichier contenant l'OTU table (TSV) de la liste d'echantillon a comparer, realise avec 'taxo_to_otutable.py', a convertir en fichier BIOM (JSON)
#2 - le fichier de metadata
	#exemple => map.tsv:
	#SampleID	Group   HealthCondition	SampleType	Date
	#sample1	1	severe	Blood 	2010-01-22
	#sample2	3	light	Blood   2010-02-16
#3 - Le fichier de sortie, avec l'extension .biom

	
#Calling example : otutable_to_biom.sh -i comp_samples_otutable.tsv -o comp_samples_otutable.biom -m map.tsv


#     .			THIS PROGRAM REQUIRES THAT BIOM HAS BEEN INSTALLED ON YOUR SESSION
#    / \		Follow the instructions below as an example of installation on bic cluster:
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

biom="$HOME/.local/bin/biom"

echo "... Converting OTU table to Biom format ..."
${biom} convert -i ${OTUTABLEFILE} -o ${OUTPUTFILE}.tmp --to-json --table-type="OTU table" --process-obs-metadata taxonomy
echo "... Adding metadata ..."
${biom} add-metadata -i ${OUTPUTFILE}.tmp -o ${OUTPUTFILE} -m ${METADATAFILE} --output-as-json
rm ${OUTPUTFILE}.tmp
