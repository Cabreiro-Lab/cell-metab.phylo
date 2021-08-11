# Phylogenetic representation

This repository contains the scripts and data to produce the tree from the paper **Increased fidelity of protein synthesis extends lifespan**, appeared in *Cell Metabolism*. 

The scripts can be run per separate, or as a whole by running the pipeline script. These scripts have been run in an Linux (Ubuntu 20.04) environment, so running them in a Windows or Mac environment might need a few tweeks. 

## Requirements

**system**

1. `mafft (>=7.460)`
2. `IQ-TREE (>=1.6.9)`

**Python**:

1. `Python >= 3.6`
2. `Biopython`

**R scripts**:
1. `optparse`
2. `seqinr`
3. `ggtree`
4. `R.devices`

### Optional!

If you want to draw a finer tree with color representation for clades and groups, there are extra requirements for the R script:

1. `ggplot2`
2. `dplyr`
3. `treeio`
4. `tidytree`

## Example of use

For the automatic pipeline, go to the scripts folder and run this command in the terminal:

```console
./pipeline.sh
```

That will make the alignments, tree and pattern detection, and will draw a simple tree within the *main_tree* folder. 

If you want to draw the other tree, run the R script in the *main_tree* folder as follows:

```console
Rscript ../scripts/drawtree_fine.R -t all_seqs.aln.treefile -p 74
```

## scripts

A simple description of each file:

1. **pipeline.sh**: main pipeline in bash
2. **pat_find.py**: Python script to detect where the aa of interest is in the alignments
3. **drawtree.R**: R script to draw the simple tree
4. **drawtree_fine.R**: R script to draw a more complete phylogenetic tree