#!/usr/bin/env Rscript

# libraries
suppressMessages(suppressWarnings(library(optparse)));
suppressMessages(library(seqinr));
suppressMessages(suppressWarnings(library(ggtree)));
library(ggplot2);
suppressMessages(library(dplyr));
suppressMessages(library(treeio));
suppressMessages(library(tidytree));
options(warn=-1)
library(R.devices)

# option list
option_list = list(
	make_option(c("-t", "--tree"), type = "character", default = NULL, 
		help = "phylogenetic tree in newik format", metavar = "character"),
	make_option(c("-p", "--position"), type = "integer", default = NULL, 
		help = "position of aminoacid", metavar = "integer")
); 

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# read the alignment
cat(paste('Reading alignment', '\n', sep = ' '))

a.aln = read.alignment('Archaea.aln', "fasta", forceToLower = TRUE)
b.aln = read.alignment('Bacteria.aln', "fasta", forceToLower = TRUE)
e.aln = read.alignment('Eukaryote.aln', "fasta", forceToLower = TRUE)
t.aln = read.alignment('all_seqs.aln', "fasta", forceToLower = TRUE)

sp_names = c(a.aln$nam, b.aln$nam, e.aln$nam)
df = data.frame(sp_names)
colnames(df) = 'label'

df = df %>% 
	mutate(Domain = ifelse(label %in% a.aln$nam, 'Archaea', 
		ifelse(label %in% b.aln$nam, 'Bacteria', 'Eukarya')))

# load trees
cat(paste('Reading tree', as.character(opt$tree), '\n',sep = ' '))

tree = read.tree(opt$tree)


### tree 1: the short tree
# simple loop to extract the AA
cat(paste('Setting AA position to', as.character(opt$position),'\n', sep = ' '))

mut = c()
for (i in 1:length(t.aln$nam)){ mut = c(mut, substr(t.aln$seq[[i]][1], 74,74)) }

# creates a dataframe with names and mutation
df2 = data.frame(t.aln$nam, mut)
colnames(df2) = c('label', 'mut')
df = df %>% left_join(df2)
df3 = df %>% rename(group = Domain) %>% as_tibble



# split names into two different groups 
groupInfo = split(df$label, df$Domain)

# generate a tree
tree1 = groupOTU(tree, groupInfo)

comp_tree = full_join(as_tibble(tree1), df3, by = "label") %>% 
select(-group.y) %>%
rename(group = group.x) %>%
as.treedata


sp = as_tibble(comp_tree)$label[1:75]

# set species names
species = c("Methanothrix thermoacetophila", "Methanococcoides burtonii", "Methanosarcina barkeri", "Methanosarcina acetivorans", "Methanosarcina mazei", 
	"Methanocorpusculum labreanum", "Methanoregula boonei", "Methanospirillum hungatei", "Methanoculleus marisnigri", "Halobacterium salinarum", 
	"Haloquadratum walsbyi", "Haloarcula marismortui", "Natronomonas pharaonis", "Halococcus morrhuae", "Archaeoglobus fulgidus", "Picrophilus torridus",
	"Thermoplasma volcanium", "Thermoplasma acidophilum", "Methanocella arvoryzae", "Schizosaccharomyces pombe", "Saccharomyces cerevisiae", "Candida glabrata", 
	"Saccharomyces uvarum", "Naumovozyma castellii", "Neurospora crassa", "Neosartorya fumigata", "Homo sapiens", "Bos taurus Bovine", "Rattus norvegicus", 
	"Chinchilla lanigera", "Sus scrofa", "Ictalurus punctatus", "Mus musculus", "Gillichthys mirabilis", "Ciona intestinalis", "Papilio dardanus", 
	"Spodoptera frugiperda", "Drosophila melanogaster", "Dermacentor variabilis", "Lumbricus rubellus", "Brugia malayi", "Caenorhabditis elegans", 
	"Tetrahymena thermophila", "Fragaria ananassa", "Euphorbia esula", "Arabidopsis thaliana", "Encephalitozoon cuniculi", "Methanobrevibacter smithii", 
	"Methanosphaera stadtmanae", "Methanothermobacter thermautotrophicus", "Pyrobaculum islandicum", "Pyrobaculum arsenaticum", "Pyrobaculum aerophilum", 
	"Pyrobaculum calidifontis", "Thermofilum pendens", "Sulfolobus acidocaldarius", "Sulfurisphaera tokodaii", "Saccharolobus solfataricus", "Hyperthermus butylicus", 
	"Aeropyrum pernix", "Staphylothermus marinus", "Ignicoccus hospitalis", "Nanoarchaeum equitans", "Methanococcus maripaludis", "Methanococcus aeolicus", 
	"Methanococcus vannielii", "Methanocaldococcus jannaschii", "Methanopyrus kandleri", "Pyrococcus horikoshii", "Pyrococcus abyssi", "Pyrococcus furiosus", 
	"Thermococcus celer", "Thermococcus kodakarensis", "Escherichia coli", "Pseudomonas putida")


new_tree = full_join(as_tibble(tree1), df3, by ="label") %>% 
select(-group.y) %>%
rename(group = group.x) 

# puting the correct species names
new_tree$label[1:75] = species

new_tree = new_tree %>%
as.treedata

# generate a tree
cat('Drawing tree... \n')
tree_plot = ggtree(new_tree,  aes(color = group) , layout = 'fan', branch.length = "none", size = 0.7) + 
	geom_tiplab(size = 1.5, aes(angle = angle, colour = mut, show.legend = FALSE)) +
	scale_colour_manual(values = c('#E8A33C',         # Archaea
								   '#37BA4B',         # Bacteria
								   '#6F7CC4',         # Eukarya
	 							   'grey60',          # K
	 							   'red'))   +        # R
	theme(legend.position = "bottom",
		legend.title = element_blank()) 

suppressGraphics(ggsave(file = 'fine_tree.pdf', tree_plot, width = 120, height = 120, units = 'mm'))

cat('Process completed! \n')