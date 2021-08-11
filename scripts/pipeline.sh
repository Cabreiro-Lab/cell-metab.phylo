#!/bin/bash
# bash script to produce a simple tree

sleep 1

printf "Running mafft to align sequences \n"

mafft --auto --quiet ../main_tree/Archaea.fasta > ../main_tree/Archaea.aln 
mafft --auto --quiet ../main_tree/Eukaryote.fasta > ../main_tree/Eukaryote.aln 
mafft --auto --quiet ../main_tree/Bacteria.fasta > ../main_tree/Bacteria.aln 
mafft --auto --quiet ../main_tree/all_seqs.fasta > ../main_tree/all_seqs.aln 

printf "Generating the tree files with IQtree"
iqtree -s ../main_tree/all_seqs.aln -st AA -nt 6 -m LG+G4 -n 0 -o Escherichia_coli -quiet

printf "Running pattern finding in the main alignment... \n"
result=$(python pat_find.py ../main_tree/all_seqs.aln)

printf "The aa is in position $result\n"

printf "Drawing simple tree!\n"
Rscript drawtree.R -a ../main_tree/all_seqs.aln -t ../main_tree/all_seqs.aln.treefile -p $result
mv simple_tree.pdf ../main_tree/

printf "Done!\n"