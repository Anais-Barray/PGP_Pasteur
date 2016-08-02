#!/bin/bash

#programme qui prend en entree :
#1 - le fichier fasta/fastq de toutes les sequences
#2 - le fichier contenant les labels du resultat de classification Kraken
#3 - la taxonomie recherchee
#4 - le fichier de sortie
#5 - l'option -r si l'utilisateur souhaite realiser l'operation inverse et recuperer les sequences qui n'ont pas la taxonomie

#note : si plusieurs taxonomies sont donnees en entree, l'ecrire comme suit : tax1\|tax2\|tax3. Dans une taxo, les espaces sont permis.

#le programme cree un fichier fasta/fastq avec uniquement les sequences de la taxonomie d'interet (ou l'inverse si specifie).

REMOVE=no

while getopts f:k:t:o:r option
do
  case $option in
  f)	FASTFILE=$OPTARG
		if [ ! -e $FASTFILE ]; then echo "   problem with input file (option -f): '$FASTFILE' does not exist" ; exit ; fi;;

  k)	KRAKENREPORT=$OPTARG
		if [ ! -e $KRAKENREPORT ]; then echo "   problem with input file (option -k): '$KRAKENREPORT' does not exist" ; exit ; fi;;
	
  t) 	TAXONOMY=$OPTARG;;

  o) 	OUTPUTFILE=$OPTARG;;

  r)	REMOVE=yes;;

  esac
done

INFILE=`basename ${FASTFILE}`;
INKRAKEN=`basename ${KRAKENREPORT}`;
DIR=`dirname ${FASTFILE}`;
fastchar=`head -c 1 $FASTFILE`;

if [ ${fastchar} = "@" ]
then
	fasttype="fastq"
	printf "\nfastq file is < $FASTFILE >"

elif [ ${fastchar} = ">" ]
then
	fasttype="fasta"
	printf "\nfasta file is < $FASTFILE >"
else
	echo file format not recognized
	exit 1
fi

printf "\nkraken label file is < $KRAKENREPORT >\n"
echo remove = $REMOVE
echo taxonomy = $TAXONOMY

echo making taxonomy lists ...
tax=`echo ${TAXONOMY} | sed "s: :_:g"` ;
TAXLIST=${DIR}/${tax}_${INKRAKEN}_list.txt ;
egrep -i "${TAXONOMY}" $KRAKENREPORT | cut -f 1 > $TAXLIST ;


if [ $REMOVE = "yes" ]	
then
	echo making notaxonomy lists ...
	NOTAXLIST=${DIR}/no${tax}_${INKRAKEN}_list.txt
	ALLTAXLIST=${DIR}/alltax_${INKRAKEN}_list.txt
	cut -f 1 $KRAKENREPORT > $ALLTAXLIST ;
	diff $TAXLIST $ALLTAXLIST | tail -n +2 | tr -d "> " > $NOTAXLIST ;
	
	if [ ${fasttype} = "fasta" ]
	then
		echo making fasta file without $TAXONOMY sequences ...
		python $SCRIPTDIR/get_fasta_from_name_seq.py $FASTFILE $NOTAXLIST $OUTPUTFILE 
	else
		echo making fastq file without $TAXONOMY sequences ...
		python $SCRIPTDIR/get_fastq_from_name_seq.py $FASTFILE $NOTAXLIST $OUTPUTFILE 
	fi
else
	if [ ${fasttype} = "fasta" ]
	then
		echo making fasta file with $TAXONOMY sequences ...
		python $SCRIPTDIR/Scripts/get_fasta_from_name_seq.py $FASTFILE $TAXLIST $OUTPUTFILE 
	else
		echo making fastq file with $TAXONOMY sequences ...
		python $SCRIPTDIR/Scripts/get_fastq_from_name_seq.py $FASTFILE $TAXLIST $OUTPUTFILE 
	fi
fi

printf "< $OUTPUTFILE > created\n"
echo erasing intermediate files ...
rm -f ${TAXLIST} ${NOTAXLIST} ${ALLTAXLIST}

echo done
