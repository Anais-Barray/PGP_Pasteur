#!/usr/bin/python3

from Bio import SeqIO
import sys, os

#options
input_file = sys.argv[1]
input_name = os.path.splitext(input_file)[0]
length_opt = int(sys.argv[2])

output_short = open(input_name+"_short.fasta", "w")
output_long = open(input_name+"_long.fasta", "w")


for record in SeqIO.parse(open(input_file), "fasta") :
	#~ print(">"+record.id)
	#~ print(record.seq)
	if (len(record.seq) < length_opt) :
		SeqIO.write(record, output_short, "fasta")
	else :
		SeqIO.write(record, output_long, "fasta")

		
output_long.close()
output_short.close()

