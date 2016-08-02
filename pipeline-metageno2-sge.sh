#!/bin/bash

### METAGENOMIC PIPELINE 2.0 : From a .bam sequence file, this pipeline will realize a cleaning, decontamination, classification and get the viral taxonomic content.
#Suited for a SGE grid parallelization
#This program has the following entries :
#1 - the input file containing the sequences to analyze (mandatory). Pipeline starts from raw .bam file, but can be relaunched from renamed .bam afterwards (debug).
#2 - an output name (mandatory)
#3 - the filtering length for Fastq Cleaner (by default = 30)
#4 - the filtering quality treshold (Phred score) for Fastq Cleaner (by default = 20)
#5 - the filtering % of bases meeting the quality threshold for Fastq Cleaner (by default = 80)

#Calling example : pipeline-metageno2-sge.sh -i sequences.bam -o Run_300516_RDC7 -l 30 -q 18 -p 70

#########################################
########## OPTIONS & VARIABLES ##########
#########################################

FQL=30 # Fqcleaner length filter
FQQ=20 # Fqcleaner Phred filter
FQP=80 # Fqcleaner % filter

usage="Usage : $(basename "$0") [-h] [-i input_file] [-o output_name] [(optional) -l -q -p length/phred/% filtering parameters]"

while getopts :i:o:l:q:p:h option
do
  case $option in
  i)	INPUT=$OPTARG
		iflag=true
		if [ ! -e $INPUT ]; then echo "err : problem with input file (option -i): '$INPUT' does not exist" ; exit ; fi;;
  o)	OUTPUT=$OPTARG
		oflag=true;;
  l)	FQL=$OPTARG;;
  q)	FQQ=$OPTARG;;
  p)	FQP=$OPTARG;;
  h)	echo "$usage"
		exit 1;;
  :)	printf "missing argument for -%s\n" "$OPTARG" >&2
		echo "$usage" >&2
		exit 1 ;;
  \?)	printf "illegal option: -%s\n" "$OPTARG" >&2
		echo "$usage" >&2
		exit 1 ;;
  esac
done
shift "$((OPTIND-1))"

if [[ -z $iflag || -z $oflag ]]
then
	echo "Exit : Missing mandatory option(s) to run the pipeline" >&2
	echo "$usage";
	exit 1;
fi

INPUT_DIR=$( cd "$(dirname "${INPUT}")" && pwd )
INPUT_NAME=$( basename "${INPUT}" )
SCRIPTDIR= :-)
KAIJUDB=${SCRIPTDIR}/kaiju/nrGenBankdb/full_parasite_db
KAIJUDIR=${SCRIPTDIR}/kaiju/bin
NTDB= :-)

#########################################
############### ANALYSIS ################
#########################################

echo "########## STARTING METAGENOMIC ANALYSIS ##########";

if [ "${INPUT_NAME}" == "${OUTPUT}.bam" ]	#if we give the renamed .bam file
then
	cd ${INPUT_DIR}
else
	if [ -d "${INPUT_DIR}/${OUTPUT}" ]	#if we give the original .bam file but the new directory exists
	then
		cd ${INPUT_DIR}/${OUTPUT}
	else 		#creates new directory and filename for the analyses
		mkdir ${INPUT_DIR}/${OUTPUT}
		cd ${INPUT_DIR}/${OUTPUT}
		ln -fs ${INPUT_DIR}/${INPUT_NAME} ${OUTPUT}.bam 
	fi
fi

# BAM 2 FASTQ
echo "... Step 1 : Bam to Fastq ..." ;
qrsh -V -cwd -now n ${SCRIPTDIR}/bam2fastq-1.1.0/bam2fastq -o ${OUTPUT}.fastq ${OUTPUT}.bam
echo
 
# FASTQCLEANER
echo "... Step 2.1 : Fastq Cleaner ..." ;
OUTPUT_CLEAN=${OUTPUT}.q${FQQ}p${FQP}l${FQL}.CFQD
qrsh -V -cwd -now n ${SCRIPTDIR}/fqCleaner.sh -f ${OUTPUT}.fastq -q ${FQQ} -l ${FQL} -p ${FQP} -x ${OUTPUT_CLEAN}.fq
echo
 
# FASTQ 2 FASTA
echo "... Step 2.2 : Fastq to Fasta ..." ;
qrsh -V -cwd -now n perl ${SCRIPTDIR}/fastq2fasta.pl ${OUTPUT_CLEAN}.fq ${OUTPUT_CLEAN}.fasta 
echo
 
# DECONSEQ
echo "... Step 3 : DeconSeq decontamination ..." ;
OUTPUT_DECONT=${OUTPUT_CLEAN}_clean
${SCRIPTDIR}/DeconSeq/pgpDeconseq_sge.sh -f ${OUTPUT_CLEAN}.fasta -d human -r viral,bacteria -c 95 -i 94 -g 1 -n 100
rm ${OUTPUT_CLEAN}_cont.fasta
echo
 
# KAIJU
echo "... Step 4.1 : Kaiju Classification ..." ;
OUTPUT_KAIJU=${OUTPUT_DECONT}.greedym1_full.kaiju
qrsh -cwd -now n -pe thread 10 -l mem_total=58G -V ${KAIJUDIR}/kaiju -i ${OUTPUT_DECONT}.fasta -o ${OUTPUT_KAIJU} -t ${KAIJUDB}/nodes.dmp -f ${KAIJUDB}/kaiju_db_nr_euk.fmi -z 10 -a greedy -m 1 -v
echo
 
echo "... Step 4.2 : Get labels ..." ;
OUTPUT_LABELS=${OUTPUT_KAIJU}.full_labels
qrsh -cwd -now n -V ${KAIJUDIR}/addTaxonNames -i ${OUTPUT_KAIJU} -o ${OUTPUT_LABELS} -t ${KAIJUDB}/nodes.dmp -n ${KAIJUDB}/names.dmp -u -p
ls *.full_labels
echo
 
# GET VIRUS FASTA SEQUENCES AND COUNTING
echo "... Step 5.1 : Get Virus fasta sequences ..." ;
OUTPUT_VIRUS=${OUTPUT_LABELS}.virus
qrsh -cwd -V -now n ${SCRIPTDIR}/get_seq_from_kaiju_taxo.sh -f ${OUTPUT_DECONT}.fasta -k ${OUTPUT_LABELS} -o ${OUTPUT_VIRUS}.fasta -t 'Virus'
echo
 
echo "... Step 5.2 : Count Virus fasta sequences ..." ;
nbseq=`grep -c "^>" ${OUTPUT_VIRUS}.fasta`;
echo "nbseq = ${nbseq}";
if [ $nbseq -gt 1000000 ]
then 
	echo "nbseq > 1 000 000, need to assemble sequences"; 
	echo
	echo "... Step 5.3 : CLC de novo assembly  ..." ;
	perl ${SCRIPTDIR}/clc4_denovo.pl -q ${OUTPUT_VIRUS}.fasta;
	cd ${OUTPUT_VIRUS}-CLC;
fi
echo

# BLASTN
echo "... Step 6.1 : Blastn on nt ..." ;
if [ $nbseq -gt 1000000 ]
then 
	OUTPUT_CONTIG=${OUTPUT_VIRUS}.contigs
	OUTPUT_UNASSEMBLED=${OUTPUT_VIRUS}.unassembled_reads
	mv contigs.fasta ${OUTPUT_CONTIG}.fasta
	mv unassembled_reads.fasta ${OUTPUT_UNASSEMBLED}.fasta
	perl ${SCRIPTDIR}/PgpBlastall/pgp_blastall -i ${OUTPUT_CONTIG}.fasta -o ${OUTPUT_CONTIG}.blastn -v 1 -b 1 -m 8 -F F -e 10E-5 -d ${NTDB} -p blastn > /dev/null 2>&1
	perl ${SCRIPTDIR}/PgpBlastall/pgp_blastall -i ${OUTPUT_UNASSEMBLED}.fasta -o ${OUTPUT_UNASSEMBLED}.blastn -v 1 -b 1 -m 8 -F F -e 10E-5 -d ${NTDB} -p blastn > /dev/null 2>&1
else
	perl ${SCRIPTDIR}/PgpBlastall/pgp_blastall -i ${OUTPUT_VIRUS}.fasta -o ${OUTPUT_VIRUS}.blastn -v 1 -b 1 -m 8 -F F -e 10E-5 -d ${NTDB} -p blastn > /dev/null 2>&1
fi
ls *.blastn
echo

echo "... Step 6.2 : Get blastn best score ..." ;
if [ $nbseq -gt 1000000 ]
then 
	OUTPUT_CONTIG_BLASTN_BEST=${OUTPUT_CONTIG}.blastn.best_score
	OUTPUT_UNASSEMBLED_BLASTN_BEST=${OUTPUT_UNASSEMBLED}.blastn.best_score
	qrsh -V -cwd -now n perl ${SCRIPTDIR}/best_score_v4.pl -i ${OUTPUT_CONTIG}.blastn -o ${OUTPUT_CONTIG_BLASTN_BEST}
	qrsh -V -cwd -now n perl ${SCRIPTDIR}/best_score_v4.pl -i ${OUTPUT_UNASSEMBLED}.blastn -o ${OUTPUT_UNASSEMBLED_BLASTN_BEST}
else
	OUTPUT_VIRUS_BLASTN_BEST=${OUTPUT_VIRUS}.blastn.best_score
	qrsh -V -cwd -now n perl ${SCRIPTDIR}/best_score_v4.pl -i ${OUTPUT_VIRUS}.blastn -o ${OUTPUT_VIRUS_BLASTN_BEST}
fi
ls *.best_score
echo

echo "... Step 6.3 : Get blastn TAXO ..." ;
if [ $nbseq -gt 1000000 ]
then 
	sh ${SCRIPTDIR}/Pf8OgsTaxoptimizer/pf8_ogs_taxoptimizer.sh ${OUTPUT_CONTIG_BLASTN_BEST}
	sh ${SCRIPTDIR}/Pf8OgsTaxoptimizer/pf8_ogs_taxoptimizer.sh ${OUTPUT_UNASSEMBLED_BLASTN_BEST}
else
	sh ${SCRIPTDIR}/Pf8OgsTaxoptimizer/pf8_ogs_taxoptimizer.sh ${OUTPUT_VIRUS_BLASTN_BEST}
fi
ls *.taxo
echo

# GET VIRUS FINAL TAXO AND CREATE ARCHIVE WITH EXCEL FILES
echo "... Step 7.1 : Creation of Virus tabulate output file(s) ..." ;

LINE="Query	Hit	%Homology_alignment	Length	MisMatch	Gap	Query_start	Query_end	Hit_start	Hit_end	E-value	Bitscore	ID	Taxo";
if [ $nbseq -gt 1000000 ]
then 
	OUTPUT_FILE_CONTIG=${OUTPUT}.contigs.blastn.virus.xls
	OUTPUT_FILE_UNASSEMBLED=${OUTPUT}.unassembled_reads.blastn.virus.xls
	LINE_CONTIG="Query	Hit	%Homology_alignment	Length	MisMatch	Gap	Query_start	Query_end	Hit_start	Hit_end	E-value	Bitscore	ID	Taxo	Nb_reads";
	
	grep -i "Virus" ${OUTPUT_CONTIG_BLASTN_BEST}.taxo > ${OUTPUT_FILE_CONTIG}.tmp
	cat ${OUTPUT_FILE_CONTIG}.tmp | cut -f 1 | sed -e 's/contig.*nbreads_//' -e 's/_.*//' > nb_read.txt
	paste -d "    " ${OUTPUT_FILE_CONTIG}.tmp nb_read.txt > ${OUTPUT_FILE_CONTIG}
	echo "${LINE_CONTIG}" > ${OUTPUT_FILE_CONTIG}.tmp
	cat ${OUTPUT_FILE_CONTIG} >> ${OUTPUT_FILE_CONTIG}.tmp
	mv ${OUTPUT_FILE_CONTIG}.tmp ${OUTPUT_FILE_CONTIG}
	
	grep -i "Virus" ${OUTPUT_UNASSEMBLED_BLASTN_BEST}.taxo > ${OUTPUT_FILE_UNASSEMBLED}
	echo "${LINE}" > ${OUTPUT_FILE_UNASSEMBLED}.tmp
	cat ${OUTPUT_FILE_UNASSEMBLED} >> ${OUTPUT_FILE_UNASSEMBLED}.tmp
	mv ${OUTPUT_FILE_UNASSEMBLED}.tmp ${OUTPUT_FILE_UNASSEMBLED}
	
	rm nb_read.txt
else
	OUTPUT_FILE=${OUTPUT}.blastn.virus.xls
	grep -i "Virus" ${OUTPUT_VIRUS_BLASTN_BEST}.taxo > ${OUTPUT_FILE}
	echo "${LINE}" > ${OUTPUT_FILE}.tmp
	cat ${OUTPUT_FILE} >> ${OUTPUT_FILE}.tmp
	mv ${OUTPUT_FILE}.tmp ${OUTPUT_FILE}
fi
ls *.xls
echo

echo "... Step 7.2 : Creation of archive ..." ;
${SCRIPTDIR}/p7zip/bin/7z a ${OUTPUT}.7z *.xls
rm *.xls
echo

echo "########## END OF METAGENOMIC ANALYSIS ##########";
