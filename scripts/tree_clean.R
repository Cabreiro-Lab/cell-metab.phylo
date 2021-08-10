## R script to read alignment and generate the phylogenetic tree
# clean version 

# libraries


library(dplyr)
library(ggplot2)
library(seqinr)
library(ggtree)
library(openxlsx)
library(treeio)
library(tidytree)
library(here)

# session options
options(width = 220)


# my library of functions
source('/Users/dmarti14/Documents/MRC_Postdoc/scripts/R_functions/all_functions.R')


# let's create a dataframe with all species within this tree

a.aln = read.alignment('Archaea_n_short.aln', "fasta", forceToLower = TRUE)
b.aln = read.alignment('Bacteria_n.aln', "fasta", forceToLower = TRUE)
e.aln = read.alignment('Eukarya_n_short.aln', "fasta", forceToLower = TRUE)
t.aln = read.alignment('tree_of_life_short.aln', "fasta", forceToLower = TRUE)

sp_names = c(a.aln$nam, b.aln$nam, e.aln$nam)
df = data.frame(sp_names)
colnames(df) = 'label'

df = df %>% 
	mutate(Domain = ifelse(label %in% a.aln$nam, 'Archaea', 
		ifelse(label %in% b.aln$nam, 'Bacteria', 'Eukarya')))


mut = c()
for (i in 1:length(t.aln$nam)){ mut = c(mut, substr(t.aln$seq[[i]][1], 172,172)) }
# creates a dataframe with names and mutation
df2 = data.frame(t.aln$nam, mut)
colnames(df2) = c('label', 'mut')

df = df %>% left_join(df2)

df3 = df %>% rename(group = Domain) %>% as_tibble



# load trees
tree = read.tree('tree_files/tree_of_life_short.aln.treefile')

root_tree = phytools::reroot(as.phylo(tree), node = 1)


# split names into two different groups 
groupInfo = split(df$label, df$Domain)

# generate a tree
tree1 = groupOTU(root_tree, groupInfo)



comp_tree = full_join(as_tibble(tree1), df3, by ="label") %>% 
select(-group.y) %>%
rename(group = group.x) %>%
as.treedata

# short version of the tree with long names
ggtree(comp_tree,  aes(color = group) , layout = 'fan', branch.length = "none", size = 0.3) + 
	geom_tiplab(size = 1.5, aes(angle = angle, colour = mut, show.legend = FALSE)) +
	# scale_color_brewer("group", palette = "Set1", type = "seq") +
	# scale_colour_manual(values = colsel(9, palette = 'pastel1')) +
	scale_colour_manual(values = c('blue',         # deletion
								   '#E8A33C',         # root
								   '#E8A33C',         # archaea
								   '#37BA4B',         # bacteria
								   'violet',         # mutation E
								   '#6F7CC4',         # eukarya
								   'grey60',         # mutation K
								   'red',         # mutation R
								   'orange')) +    # mutation T
	theme(legend.position = "bottom",
		legend.title = element_blank()) 

ggsave(file = here('analysis', 'short_tree_v1.pdf'), width = 150, height = 150, units = 'mm', scale = 2, device = 'pdf')




# get the names of the species

sp = as_tibble(comp_tree)$label[1:501]


species = c("Escherichia coli", "Pseudomonas putida", "Natrialba hulunbeirensis", "Haloterrigena daqingensis", "Natronorubrum tibetense", "Natrialba asiatica st ATCC 700177", "Natrinema gari", 
  "Haloterrigena turkmenica st ATCC 51198", "Natronolimnobius baerhuensis", "Natronorubrum texcoconense", "Natronolimnobius sp AArcMg", "Natronolimnobius aegyptiacus", 
  "Natronococcus amylolyticus", "Halococcus hamelinensis", "Halarchaeum acidiphilum", "Halobacteriales archaeon SW_5_68_122", "Haloarcula rubripromontorii", 
  "Haloarcula amylolytica JCM 13557", "Halobacteriales archaeon QH_8_67_27", "Haloarcula argentinensis", "Haladaptatus litoreus", "Halobacteriales archaeon QS_1_67_19", 
  "Halorientalis sp IM1011", "Haloferax mucosum", "Haloferax sp Atlit6N", "Haloferax lucentense strain DSM 14919", "Haloferax sp Atlit48N", "Haloferax sulfurifontis ATCC", 
  "Halorubrum lacusprofundi strain ATCC 49239", "Halorubrum sp 48_1_W", "Halorubrum sp C191", "Halorubrum tropicale", "Halorubrum distributum", "Halorubrum vacuolatum",
  "Halopenitus malekzadehii", "Haloferax sp Atlit24N", "Halobacterium jilantaiense", "Haloquadratum walsbyi strain DSM 16790 ", "Halobacteriales archaeon QH_6_64_20", 
  "Halobacteriales archaeon QH_8_64_26", "Halobacteriales archaeon SW_12_71_31", "Halobaculum gomorrense", "Halobacteriales archaeon SW_8_66_22", "Halococcus saccharolyticus DSM 5350", 
  "Halobacteriales archaeon SW_9_67_24", "Halobacteriales archaeon SW_7_65_23", "Halovivax ruber strain DSM 18193", "Halobacteriales archaeon SW_12_69_24", 
  "Halobacteriales archaeon QH_3_68_24", "Nanohaloarchaea archaeon SW_7_43_1", "Nanosalina sp strain J07AB43", "Nanohaloarchaea archaeon SG9", 
  "Nanohaloarchaea archaeon SW_7_46_7", "Haloredivivus sp strain G17 H0AEB7", "archaeon SCGAAA382B04", "Methanosarcina barkeri 227", "Methanosarcina sp Ant1", 
  "Methanosarcina mazei SarPi", "Methanosarcina mazei LYC", "Methanosarcina sp 2HT1A15", "Methanosarcina sp 2HT1A6", "Methanosarcina sp WH1", "Methanosarcina sp 795",
   "Methanosarcina sp DSM 11855", "Methanomethylovorans sp PtaU1Bin073", "Methanomethylovorans sp PtaU1Bin093", "Methanohalophilus portucalensis FDF1", "Methanosalsum zhilinae strain DSM 4017", 
   "Methanohalobium evestigatum strain ATCC BAA1072", "Candidatus Methanoperedens nitroreducens", "Candidatus Methanoperedens sp BLZ1", "Candidatus Methanoperedenaceae archaeon 1", 
   "ANME2 cluster archaeon", "Methanoculleus bourgensis strain ATCC 43281", "Methanoculleus thermophilus", "Methanoregula sp PtaU1Bin051", "Methanoregula formicica strain DSM 22288",
  "Methanospirillum hungatei JF1 strain ATCC 27890", "Methanofollis liminatans", "Marine Group II euryarchaeote 1", "Marine Group II euryarchaeote 2", "Euryarchaeota archaeon 1", 
  "Euryarchaeota archaeon 2", "Euryarchaeota archaeon 3", "uncultured marine group IIIII SAT1000", "Euryarchaeota archaeon 4", "Euryarchaeota archaeon 5", "Marine Group II euryarchaeote 3", 
  "Euryarchaeota archaeon 6", "Euryarchaeota archaeon 7", "Euryarchaeota archaeon 8", "uncultured marine group IIIII KM3", "Marine Group II euryarchaeote 4", "Marine Group II euryarchaeote 5", 
  "Euryarchaeota archaeon 9", "Thermoplasmata archaeon", "Euryarchaeota archaeon 10", "Euryarchaeota archaeon 11", "Marine Group II euryarchaeote 6", "Acidiplasma cupricumulans", 
  "Acidiplasma sp MBA1", "Acidiplasma aeolicum", "Picrophilus torridus strain ATCC 700027", "Candidatus Methanoplasma termitum", "Methanomassiliicoccales archaeon RumEn", 
  "Euryarchaeota archaeon 12 RBG_16_67_27", "Marine Group III euryarchaeote CGEpi2", "Marine Group III euryarchaeote CGEpi1", "Marine Group III euryarchaeote CGEpi3", 
  "Euryarchaeota archaeon 13", "Euryarchaeota archaeon 14", "Euryarchaeota archaeon 15", "Candidatus Altiarchaeales archaeon HGWA1", "Candidatus Altiarchaeum sp CG_4_9_14_0_8", 
  "Candidatus Altiarchaeum sp CG_4_10_14_0_8", "Thermoplasmatales archaeon SM150", "Thermoplasmata archaeon M11B2D", "Thermoplasmatales archaeon SG8524", 
  "Archaeoglobus veneficus strain DSM 11195", "Archaeoglobus fulgidus", "ANME1 cluster archaeon", "uncultured archaeon D1JEK6", "Methanocella arvoryzae strain DSM 22066", 
  "Methanocella conradii strain DSM 24694", "Methanocella paludicola strain DSM 17711", "Ignicoccus hospitalis strain KIN4I", "Metallosphaera yellowstonensis MK1", 
  "Acidianus brierleyi", "Acidianus hospitalis strain W1", "Sulfolobus acidocaldarius", "Sulfolobus sp SCGC AB777 L09", "Sulfolobus sp SCGC AB777 G06", "Sulfolobus islandicus strain YG5714 ", 
  "Sulfolobus islandicus strain M1425 ", "Sulfolobus islandicus strain LS215 ", "archaeon HR02", "archaeon HR01", "Thaumarchaeota archaeon ex4484", "Candidatus Marsarchaeota G1 archaeon BE_D", 
  "Candidatus Marsarchaeota G1 archaeon OSP_B", "Candidatus Marsarchaeota G2 archaeon OSP_D", "Candidatus Marsarchaeota G2 archaeon ECH_B", "uncultured Acidilobus sp OSP8", 
  "Acidilobus sp SCGC AC742 E15", "Staphylothermus hellenicus strain DSM 12710", "Desulfurococcales archaeon ex4484 58", "Thermoproteus tenax strain ATCC 35583", 
  "Thermoproteus sp CP80", "Thermoproteus sp CIS19", "Thermocladium sp ECHB", "Pyrobaculum aerophilum", "Pyrobaculum sp WP30", "Pyrobaculum neutrophilum strain DSM 2338", 
  "Pyrobaculum oguniense strain DSM 13380", "Pyrobaculum arsenaticum strain DSM 13514", "Pyrobaculum ferrireducens", "Pyrobaculum aerophilum strain ATCC 51768", 
  "Candidatus Korarchaeota archaeon", "Candidatus Nanoclepta minutus", "miscellaneous Crenarchaeota group15", "Pyrococcus horikoshii strain ATCC 700860", "Pyrococcus abyssi strain GE5", 
  "Pyrococcus furiosus strain ATCC 43587", "Pyrococcus furiosus COM1", "Pyrococcus sp ST04", "Thermococcus sibiricus", "Thermococcus sibiricus strain MM 739", "Thermococcus sp 2319x1", 
  "Thermococcales archaeon 44 46", "Thermococcus pacificus", "Thermococcus cleftensis", "Methanocaldococcus jannaschii strain ATCC 43067", "Methanocaldococcus villosus KIN24T80", 
  "Methanotorris formicicus McS70", "Methanococcus aeolicus strain ATCC BAA1280", "Methanococcus maripaludis OS7", "Methanococcus maripaludis strain C5 ATCC BAA1333", "archaeon HR05", 
  "Thermofilum pendens strain DSM 2475", "Methanosphaera sp rholeuAM6", "Methanosphaera stadtmanae", "Methanobacterium sp MZA1", "Methanobacterium formicicum", "Methanobrevibacter gottschalkii", 
  "Methanobrevibacter sp YE315", "Methanobrevibacter millerae", "Methanobrevibacter curvatus", "Methanobrevibacter smithii DSM 2374", "Methanobacterium lacus strain AL21", 
  "Methanothermobacter tenebrarum", "candidate divison MSBL1 archaeon", "Hadesarchaea archaeon YNP", "Hadesarchaea archaeon DG331", "Theionarchaea archaeon DG701", 
  "Arc I group archaeon B15fssc0709 Meth", "Candidatus Aenigmarchaeota archaeon CG01", "Candidatus Aenigmarchaeota archaeon CG_410", "Candidatus Micrarchaeota archaeon CG1_02_60", 
  "Candidatus Micrarchaeota archaeon CG4_10_14", "archaeon A0A2D6WUZ7", "archaeon Candidatus Huberarchaea CG18_big_fil", "archaeon CG07 land 8_20_14", 
  "Candidatus Diapherotrites archaeon 1", "Candidatus Pacearchaeota archaeon RBG_13_33", "Archaeon GW2011 AR13", "Candidatus Diapherotrites archaeon ADurbBin253", 
  "Candidatus Pacearchaeota archaeon CG10_big_fil 82114010319", "Candidatus Pacearchaeota archaeon A0A2D6TZP8", "Candidatus Pacearchaeota archaeon C10_23914", 
  "Candidatus Pacearchaeota archaeon CG4_101402", "Candidatus Pacearchaeota archaeon CG4_101408", "Candidatus Pacearchaeota archaeon CG10_big_fil 8211401035219", 
  "Candidatus Pacearchaeota archaeon CG10_big_fil 821140103513", "Candidatus Pacearchaeota archaeon CG49143317", "Candidatus Woesearchaeota archaeon CG0882014020479", 
  "archaeon D22 A0A202DF04", "Candidatus Woesearchaeota archaeon", "Candidatus Woesearchaeota archaeon CG10 big_fil 82114013312", "Candidatus Woesearchaeota archaeon CG4101408475", 
  "Candidatus Pacearchaeota archaeon ex4484 26", "Candidatus Heimdallarchaeota archaeon LC 3", "Candidatus Thorarchaeota archaeon SMTZ45", "uncultured marine thaumarchaeote KM390G11", 
  "Candidatus Bathyarchaeota archaeon CG07 82014080479", "Candidatus Bathyarchaeota archaeon ex4484 231", "Candidatus Bathyarchaeota archaeon B261", "miscellaneous Crenarchaeota group archaeon SMZ155", 
  "Lokiarchaeum sp GC14 75", "archaeon 13120CM2549", "uncultured marine thaumarchaeote KM3 79 D07", "uncultured marine crenarchaeote HF4000 APKG5C13", "Thaumarchaeota archaeon SCGC AC337 F14", 
  "uncultured marine crenarchaeote HF4000 APKG3E18", "uncultured marine thaumarchaeote KM3 35 E05", "Candidatus Nitrosopelagicus brevis", "uncultured marine thaumarchaeote KM3 144 G01", "Candidatus Nitrosopumilus sp NM25", 
  "Candidatus Nitrosopumilus sediminis", "Marine Group I thaumarchaeote SCGC RSA3", "Candidatus Nitrosopumilus salaria BD31", "Candidatus Nitrosopumilus adriaticus", "Nitrosopumilus sp 1", "Nitrosopumilus sp 2", 
  "Candidatus Nitrosoarchaeum sp", "Candidatus Nitrosoarchaeum limnia SFB1", "Candidatus Nitrosoarchaeum limnia BG20", "Nitrosarchaeum koreense MY1", "Candidatus Nitrososphaera sp 13", 
  "Candidatus Nitrocosmicus oleophilus", "Nicotiana attenuata", "Citrus clementina", "Aquilegia coerulea", "Gossypium barbadense", "Gossypium raimondii", "Cephalotus follicularis", 
  "Dorcoceras hygrometricum", "Cuscuta australis", "Juglans regia", "Eucalyptus grandis", "Solanum chacoense", "Coffea canephora", "Lupinus angustifolius", "Nicotiana sylvestris",
   "Populus trichocarpa", "Lotus japonicus", "Cicer arietinum", "Glycine max", "Anthurium amnicola", "Triticum aestivum", "Oryza sativa subsp japonica", "Hordeum vulgare subsp vulgare", "Oryza sativa subsp indica", 
   "Oryza nivara", "Oryza punctata", "Panicum hallii var hallii", "Ananas comosus", "Zea mays", "Wollemia nobilis", "Apostasia shenzhenica", "Brassica rapa subsp pekinensis", "Arabidopsis thaliana", 
   "Corchorus capsularis Jute", "Chlamydomonas reinhardtii", "Noccaea caerulescens", "Auxenochlorella protothecoides 1", "Auxenochlorella protothecoides 2", "Ostreococcus tauri", "Hirondellea gigas", 
   "Karenia brevis", "Corethrella appendiculata", "Anopheles christyi", "Aedes albopictus", "Anopheles gambiae", "Heliconius melpomene", "Papilio polytes", "Rhodnius neglectus", "Dermacentor variabilis", 
   "Pelinobius muticus", "Tropilaelaps mercedesae", "Maconellicoccus hirsutus", "Drosophila yakuba", "Glossina palpalis gambiensis", "Haematobia irritans", "Glossina morsitans morsitans", "Drosophila erecta", 
   "Apis cerana", "Habropoda laboriosa", "Trachymyrmex zeteki", "Trichomalopsis sarcophagae", "Agrilus planipennis", "Orchesella cincta", "Capitella teleta", "Periplaneta americana", "Cuerna arida", 
   "Crassostrea gigas", "Diploscapter pachys", "Heligmosomoides polygyrus bakeri", "Toxocara canis", "Trichinella patagoniensis", "Trichinella pseudospiralis", "Macrostomum lignano", "Strongyloides ratti", 
   "Rhabditophanes sp KR3021", "Echinostoma caproni", "Schistosoma margrebowiei", "Schistosoma mansoni", "Biomphalaria glabrata", "Arion vulgaris", "Pomacea canaliculata", "Ictidomys tridecemlineatus", 
   "Papio anubis", "Cebus capucinus imitator", "Rhinopithecus roxellana", "Leptonychotes weddellii", "Propithecus coquereli", "Pongo abelii", "Lipotes vexillifer", "Tarsius syrichta", "Rhinopithecus bieti 1", 
   "Ailuropoda melanoleuca", "Rhinopithecus bieti 2", "Macaca fascicularis", "Daphnia magna", "Rhinopithecus bieti 3", "Saimiri boliviensis boliviensis", "Aotus nancymaae", "Rhinopithecus bieti 4", "Macaca mulatta", 
   "Colobus angolensis", "Callithrix jacchus", "Pan paniscus", "Pan troglodytes", "Limosa lapponica", "Rhinopithecus bieti 5", "Mus musculus", "Antrostomus carolinensis", "Sus scrofa", "Lonchura striata domestica", 
   "Solea senegalensis", "Poecilia latipinna sailfin", "Nothobranchius korthausae", "Lithobates catesbeiana", "Iconisemion striatum", "Nothobranchius kuhntae", "Stegastes partitus", "Nothobranchius furzeri", 
   "Ictalurus punctatus", "Macaca nemestrina", "Caligus rogercresseyi", "Micrurus surinamensis", "Stichopus japonicus", "Plasmodium chabaudi", "Plasmodium yoelii", "Plasmodium sp gorilla clade G2", "Plasmodium malariae", 
   "Toxoplasma gondii", "Toxoplasma gondii strain ATCC 50611", "Cryptosporidium ubiquitum", "Spizellomyces punctatus DAOM BR117", "Trichoderma gamsii", "Hypocrea jecorina strain ATCC 56765", "Trichoderma harzianum CBS 22695", 
   "Trichoderma arundinaceum", "Hypocrea jecorina strain QM6a", "Fusarium oxysporum f sp cubense strain race 1", "Gibberella zeae strain PH1  ATCC MYA4620", "Fusarium proliferatum strain ET1", 
   "Fusarium oxysporum f sp pisi", "Lomentospora prolificans", "Trichoderma parareesei", "Ophiocordyceps polyrhachisfurcata", "Tolypocladium paradoxum", "Cordyceps fumosorosea ARSEF", "Cordyceps confragosa", 
   "Stachybotrys elegans", "Trichoderma longibrachiatum ATCC 18648", "Verticillium dahliae Verticillium wilt", "Verticillium dahliae strain VdLs17", "Colletotrichum tofieldiae", "Colletotrichum orchidophilum", 
   "Sordaria macrospora strain ATCC MYA33", "Madurella mycetomatis", "Thielavia terrestris strain ATCC 38088", "Chaetomium thermophilum strain DSM 1495", "Valsa mali var pyri", "Togninia minima strain UCRPA7", 
   "Magnaporthe grisea", "Diaporthe ampelina", "Aspergillus turcosus", "Aspergillus carbonarius strain", "Aspergillus luchuensis", "Aspergillus oryzae strain ATCC 42149", "Aspergillus flavus", 
   "Aspergillus niger strain CBS 51388", "Aspergillus calidoustus", "Penicillium digitatum strain Pd1 CECT 20795", "Penicillium coprophilum", "Penicillium nordicum", "Penicilliopsis zonata CBS 50665", 
   "Talaromyces stipitatus strain ATCC 10500 ", "Talaromyces islandicus", "Rasamsonia emersonii CBS 39364", "Baudoinia panamericana strain UAMH 10762", "Zymoseptoria tritici ST99CH", 
   "Schizosaccharomyces octosporus strain yFS286", "Endocarpon pusillum strain Z07020  HMASL300199", "Emmonsia crescens UAMH", "Emmonsia sp CAC2015a", "Ajellomyces dermatitidis strain ER3", 
   "Blastomyces gilchristii strain SLH14081", "Trichophyton verrucosum strain HKI 0517", "Trichophyton rubrum", "Trichophyton rubrum strain ATCC MYA4607", "Fonsecaea multimorphosa", "Araucaria cunninghamii", 
   "Exophiala mesophila", "Cladophialophora carrionii", "Diplodia corticola", "Alternaria alternata", "Corynespora cassiicola", "Rutstroemia sp NJR2017a", "Phialocephala scopiformis", "Meliniomyces bicolor E", 
   "Phialophora cf hyalina BP 5553", "Pseudogymnoascus sp VKM", "Pseudogymnoascus destructans strain ATCC MYA485", "Erysiphe pulchra", "Glarea lozoyensis strain ATCC 20868 ", "Cyberlindnera jadinii strain ATCC 18201", 
   "Lachancea thermotolerans strain ATCC 56472", "Lachancea lanzarotensis", "Hanseniaspora guilliermondii", "Saccharomyces cerevisiae strain RM111a", "Saccharomyces cerevisiae strain Kyokai no 7", 
   "Zygosaccharomyces bailii strain CLIB 213", "Zygosaccharomyces rouxii", "Zygosaccharomyces parabailii", "Metschnikowia bicuspidata var bicuspidata NRRL", "Debaryomyces hansenii strain ATCC 36239", 
   "Candida maltosa strain Xu316 Yeast M3JT75", "Meyerozyma guilliermondii strain ATCC 6260", "Pichia sorbitophila strain ATCC MYA4447", "Pichia membranifaciens NRRL Y2026", "Pichia kudriavzevii Yeast", 
   "Lipomyces starkeyi NRRL Y11557", "Wickerhamiella sorbophila", "Medicago truncatula", "Tuber aestivum summer truffle", "Uncinocarpus reesii strain UAMH 1704", "Ustilaginoidea virens", 
   "Pochonia chlamydosporia 170", "Torrubiella hemipterigena", "Colletotrichum chlorophyti", "Sphaceloma murrayae", "Pestalotiopsis fici strain W1061 ", "Microbotryum intermedium", "Lentinula edodes", 
   "Pycnoporus cinnabarinus", "Laetiporus sulphureus 9353", "Hypholoma sublateritium FD334", "Amanita muscaria Koide BX008", "Coprinopsis cinerea strain Okayama7  130", "Malassezia globosa strain ATCC MYA4612", 
   "Fopius arisanus", "Kwoniella dejecticola CBS 10117", "Cryptococcus depauperatus CBS 7841", "Kockovaella imperatae", "Cryptococcus neoformans var grubii serotype A strain H99", "Sistotremastrum niveocremeum HHB9708", 
   "Rhizoctonia solani 123E", "Puccinia striiformis", "Puccinia triticina isolate 11 race 1 BBBD", "Melampsora laricipopulina strain 98AG31", "Basidiobolus meristosporus CBS", "Piromyces finnis", "Smittium simulii", 
   "Absidia repens", "Rhizopus microsporus", "Rhizophagus irregularis strain DAOM 181602", "Batrachochytrium dendrobatidis", "Nematocida sp 1 ERTm6")






new_tree = full_join(as_tibble(tree1), df3, by ="label") %>% 
select(-group.y) %>%
rename(group = group.x) 

# puting the correct species names
new_tree$label[1:501] = species

new_tree = new_tree %>%
as.treedata



ggtree(new_tree,  aes(color = group) , layout = 'fan', branch.length = "none", size = 0.3) + 
	geom_tiplab(size = 1.5, aes(angle = angle, colour = mut, show.legend = FALSE)) +
	# scale_color_brewer("group", palette = "Set1", type = "seq") +
	# scale_colour_manual(values = colsel(9, palette = 'pastel1')) +
	scale_colour_manual(values = c('blue',         # deletion
								   '#E8A33C',         # root
								   '#E8A33C',         # archaea
								   '#37BA4B',         # bacteria
								   'violet',         # mutation E
								   '#6F7CC4',         # eukarya
								   'grey60',         # mutation K
								   'red',         # mutation R
								   'orange')) +    # mutation T
	theme(legend.position = "bottom",
		legend.title = element_blank()) 


ggsave(file = here('analysis', 'short_tree_v2.pdf'), width = 150, height = 150, units = 'mm', scale = 2, device = 'pdf')




###################
### large trees ###
###################


# let's create a dataframe with all species within this tree

a.aln = read.alignment('Archaea_n.aln', "fasta", forceToLower = TRUE)
b.aln = read.alignment('Bacteria_n.aln', "fasta", forceToLower = TRUE)
e.aln = read.alignment('Eukarya_n.aln', "fasta", forceToLower = TRUE)
t.aln = read.alignment('tree_of_life.aln', "fasta", forceToLower = TRUE)

sp_names = c(a.aln$nam, b.aln$nam, e.aln$nam)
df = data.frame(sp_names)
colnames(df) = 'IDs'

df = df %>% 
	mutate(Domain = ifelse(IDs %in% a.aln$nam, 'Archaea', 
		ifelse(IDs %in% b.aln$nam, 'Bacteria', 'Eukarya')))


mut = c()
for (i in 1:length(t.aln$nam)){ mut = c(mut, substr(t.aln$seq[[i]][1], 539,539)) }
# creates a dataframe with names and mutation
df2 = data.frame(t.aln$nam, mut)
colnames(df2) = c('IDs', 'mut')

df = df %>% left_join(df2)

df3 = df %>% rename(label = IDs, group = Domain) %>% as_tibble



# load trees
tree = read.tree('tree_files/tree_of_life.aln.treefile')


# split names into two different groups 
groupInfo = split(df$IDs, df$Domain)

# generate a tree
tree1 = groupOTU(tree, groupInfo)


new_tree = full_join(as_tibble(tree1), df3, by ="label") %>% 
select(-group.y) %>%
rename(group = group.x) %>%
as.treedata

ggtree(new_tree,  aes(color = group) , layout = 'fan', branch.length = "none", size = 0.3) + 
	geom_tiplab(size = 0.5, aes(angle = angle, colour = mut, show.legend = FALSE)) +
	# scale_color_brewer("group", palette = "Set1", type = "seq") +
	# scale_colour_manual(values = colsel(9, palette = 'sat1')) +
	scale_colour_manual(values = c('blue',         # deletion
								   '#E8A33C',         # root
								   '#E8A33C',         # archaea
								   '#37BA4B',         # bacteria
								   'violet',         # mutation E
								   '#6F7CC4',         # eukarya
								   'grey60',         # mutation K
								   'red',         # mutation R
								   'orange')) +    # mutation T
	theme(legend.position = "bottom",
		legend.title = element_blank()) 

ggsave(file = 'large_tree.pdf', width = 250, height = 250, units = 'mm', scale = 2, device = 'pdf')



























