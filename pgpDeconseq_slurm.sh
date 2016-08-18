#!/bin/sh

#programme qui prend en entree :
#1 - le fichier fasta/fastq de toutes les sequences
#2 - la database servant a decontaminer (ex : human)
#3 - la database servant a garder les faux positifs (ex : parasitic or viral or bacteria)
#4 - le % min de couverture d'un read sur le genome de reference pour le considerer comme mappe (ex : 95)
#5 - le % min d'identité d'un read sur le genome de reference pour le considerer comme mappe (ex : 94)
#6 - la procedure de regroupement des sequences. group = 1 pour regrouper les reads cleanes et ambigus, group = 2 pour regrouper les reads contamines et ambigus.

#le programme genere deux fichiers fasta selon que les reads initiaux aient mappe ou non sur les genomes de references, et selon la facon dont ils ont ete regroupes
# usage : ./pgpDeconseq_slurm.sh -f RunA.q18p70.CFQD.fasta -d human -r parasitic -c 95 -i 94 -g 1 -n 200

### normal submission : srun -q cibu -cwd -now n -V ../deconseq.pl -f RunA.q18p70.CFQD.fasta -dbs human -dbs_retain viral,bacteria -c 95 -i 94 -group 1


BIN_DIR=/pasteur/projets/specific/PGP_Work/Scripts/DeconSeq
CWD=$PWD

### OPTIONS
while getopts f:d:r:c:i:g:n: option
do
  case $option in
  f)	FASTFILE=$OPTARG;;
  d)	DB_CONT=$OPTARG;;
  r)	DB_RETAIN=$OPTARG;;
  c)	COVERAGE=$OPTARG;;
  i)	IDENTITY=$OPTARG;;
  g)	GROUP=$OPTARG;;
  n)	NB_FILE=$OPTARG;;
  esac
done

infile=`basename ${FASTFILE}`;
indir=`dirname ${FASTFILE}`;
input_name=${infile%.*}

WORKING_DIR=${indir}/${infile}-DECONSEQ

if [[ ! -d ${WORKING_DIR} ]]; then
    echo creating ${WORKING_DIR}/;
    mkdir ${WORKING_DIR};
    cd ${WORKING_DIR}/ ;
else
    echo ${WORKING_DIR} exists! ;
    exit 1;
fi

ln -fs ../${infile}

echo splitting input-file up in ${NB_FILE} files for parallelized decontamination ...
srun ${BIN_DIR}/splitFasta.pl -verbose -i ${infile} -n ${NB_FILE} ;
echo done

echo creating sbatch submission file ...
OUTFILE=generated.sbatch.${input_name}         # Name of the file to generate.

(
cat <<EOF
#!/bin/sh

		srun ${BIN_DIR}/deconseq.pl -f ${infile}_c\${SLURM_ARRAY_TASK_ID}.fasta -dbs ${DB_CONT} -dbs_retain ${DB_RETAIN} -c ${COVERAGE} -i ${IDENTITY} -group ${GROUP} -id \${SLURM_ARRAY_JOB_ID}.\${SLURM_ARRAY_TASK_ID} ;

    exit 0
EOF
) > $OUTFILE

echo done.

echo submitting ...
THIS_JOB=`sbatch --array=1-${NB_FILE} ${OUTFILE}` ;
JOB_ID=${THIS_JOB##* } ;
echo job id is ${JOB_ID} ;
echo done submitting

echo
echo checking if job finished ...
sleep 30;
while : ; do
    count_nb_running_tasks=`squeue | grep ${JOB_ID} | wc -l`;
    [[ ${count_nb_running_tasks} -eq 0 ]] && break; 
    sleep 30;
done

echo job finished
echo
echo cleaning up ...
echo

cat *_cont.fa > ${input_name}_cont.fasta ;
cat *_clean.fa > ${input_name}_clean.fasta ;
mv ${input_name}_clean.fasta .. ;
mv ${input_name}_cont.fasta .. ;

echo
cd $CWD ;
if [[ -d ${WORKING_DIR} ]]; then
    rm -f -r ${WORKING_DIR}/ ;
else
    echo ${WORKING_DIR} is not there ;
fi

echo
echo done.
echo
echo exiting.
echo
exit 0 ;
