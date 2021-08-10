#!/usr/bin/env Rscript

# libraries
suppressMessages(suppressWarnings(library(optparse)));
suppressMessages(library(seqinr));
suppressMessages(suppressWarnings(library(ggtree)));
options(warn=-1)
suppressMessages(library(R.devices))

# option list
option_list = list(
	make_option(c("-a", "--aln"), type = "character", default = NULL, 
		help = "alignment file, in fasta format", metavar = "character"),
	make_option(c("-t", "--tree"), type = "character", default = NULL, 
		help = "phylogenetic tree in newik format", metavar = "character"),
	make_option(c("-p", "--position"), type = "integer", default = NULL, 
		help = "position of aminoacid", metavar = "integer")
); 

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# read the alignment
cat(paste('Reading alignment', as.character(opt$aln), '\n', sep = ' '))

aln = read.alignment(opt$aln, "fasta", forceToLower = TRUE)


# load trees
cat(paste('Reading tree', as.character(opt$tree), '\n',sep = ' '))

tree = read.tree(opt$tree)


### tree 1: the short tree
# simple loop to extract the AA
cat(paste('Setting AA position to', as.character(opt$position),'\n', sep = ' '))

mut = c()
for (i in 1:length(aln$nam)){ mut = c(mut, substr(aln$seq[[i]][1], opt$position, opt$position)) }

# creates a dataframe with names and mutation
df = data.frame(aln$nam, mut)

# split names into two different groups 
groupInfo = split(df$aln.nam, df$mut)



# generate a tree
cat('Drawing tree... \n')
tree1 = groupOTU(tree, groupInfo)
tree_plot = ggtree(tree1, aes(color = group), layout = 'circular', branch.length = "none") + 
	geom_tiplab(size = 2, aes(angle = angle)) + 
	theme(legend.position = c(0.54,0.455)) 

suppressGraphics(ggsave(file = 'simple_tree.pdf', tree_plot, width = 120, height = 120, units = 'mm'))

cat('Process completed! \n')