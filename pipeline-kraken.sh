#!/bin/bash

### PIPELINE DES ANALYSES KRAKEN (kraken, labels, rapport)
#programme qui prend en entree :
#1 - le fichier fasta/fastq de toutes les sequences
#2 - la base de donnees interrogee : parasitic, viral ou full
#3 - le filtrage (par defaut 0, cf Kraken manual).
#4 - le nombre de threads utilises par Kraken.

#Exemple d'appel : sh pipeline-kraken.sh -i sequences.fasta -d parasitic -f 0.51 -t 8

# For qsub purposes :
#~ #$ -S /bin/bash
#~ #$ -q cibu
#~ #$ -cwd 
#~ #$ -pe thread 8
#~ #$ -now n
#~ #$ -l mem_total=600G

### OPTIONS
while getopts i:d:f:t: option
do
  case $option in
  i)	input_file=$OPTARG
		if [ ! -e $input_file ]; then echo "err : problem with input file (option -i): '$input_file' does not exist" ; exit ; fi;;

  d)	kraken_db=$OPTARG
		if [ ! -e $SCRIPTDIR/Kraken/BUILD/databases/${kraken_db} ]; then echo "err : problem with krakenq database (option -d): '$kraken_db' does not exist" ; exit ; fi;;
		
  f)	filter=$OPTARG;;
  t)	thread=$OPTARG;;
  esac
done


### VARIABLES
input_name=${input_file%.*}
input_dir=$( dirname "$input_file" )
dir_kraken=$SCRIPTDIR/Kraken/BUILD
db=${dir_kraken}/databases/${kraken_db}
output_name=${input_dir}/${input_name}.${kraken_db}.kraken
if [ ${filter} ]; then output_name_filter=${input_dir}/${input_name}.${kraken_db}.filter${filter}.kraken; fi
if [ -z $thread ]; then echo "    no thread number specified (option -t), using 1 thread by default " ; thread=1 ; fi
qrsht="qrsh -V -q cibu -cwd -pe thread ${thread} -l mem_total=600G -now n"
qrshn="qrsh -V -q cibu -cwd -now n"

#Test file format
fastchar=`head -c 1 $input_file`
if [ ${fastchar} = "@" ]
then
	fasttype="fastq"
elif [ ${fastchar} = ">" ]
then
	fasttype="fasta"
else
	echo file format not recognized
	exit 1
fi


### KRAKEN
#Classification
if [ ! -e $output_name ]
then
	echo making Kraken classification ...
	if [ ${fasttype} = "fastq" ] 
	then 
		${qrsht} ${dir_kraken}/kraken --threads $thread --db $db --fastq-input ${input_file} > ${output_name}
	else
		${qrsht} ${dir_kraken}/kraken --threads $thread --db $db ${input_file} > ${output_name}
	fi
fi

#Filter
if [ $filter ]
then
	echo filtering Kraken matches ...
	${qrshn} ${dir_kraken}/kraken-filter --db $db --threshold $filter ${output_name} > ${output_name_filter}
	echo assigning labels \(NCBI taxonomy\) to matches ...
	${qrshn} ${dir_kraken}/kraken-translate --db $db ${output_name_filter} > ${output_name_filter}.labels 	#Labels
	echo creating Kraken report ...
	${qrshn} ${dir_kraken}/kraken-report --db $db ${output_name_filter} > ${output_name_filter}.report		#Report
else
	echo assigning labels \(NCBI taxonomy\) to matches ...
	${qrshn} ${dir_kraken}/kraken-translate --db $db ${output_name} > ${output_name}.labels 	#Labels
	echo creating Kraken report ...
	${qrshn} ${dir_kraken}/kraken-report --db $db ${output_name} > ${output_name}.report 	#Report
fi

