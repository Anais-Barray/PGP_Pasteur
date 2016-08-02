#!/usr/bin/python3

#programme qui prend en entree :
#1 - le fichier fastq de toutes les sequences
#2 - le fichier contenant les ID des sequences a recuperer
#3 - le fichier de sortie

#le programme cree un fichier fastq avec uniquement les sequences d'interet

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import sys

input_file = sys.argv[1]
id_file = sys.argv[2]
output_file = sys.argv[3]

wanted = set(line.rstrip("\n").split(None,1)[0] for line in open(id_file))
print "Found %i unique identifiers in %s" % (len(wanted), id_file)

count = 0
handle = open(output_file, "w")

for title, seq, qual in FastqGeneralIterator(open(input_file)) :
    if title.split(None,1)[0] in wanted:
        handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        count += 1
handle.close()

print "Saved %i records from %s to %s" % (count, input_file, output_file)
if count < len(wanted):
    print "Warning %i IDs not found in %s" % (len(wanted)-count, input_file)
