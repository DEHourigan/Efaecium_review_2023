library(ape)
library(msa)
library(ape)
library(ggtree)
library(seqinr)
library(ggplot2)


# Load your protein sequences from the FASTA file
ligases <- readAAStringSet("/home/david/ligases/ligases_labelled.faa")

# Perform multiple sequence alignment using MUSCLE algorithm
ligases_msa <- msa(ligases, method="Muscle" )
msaPrettyPrint(ligases_msa, file="/home/david/ligases_labelled.tex", output="tex",
 alFile = "/home/david/ligases/alignment.fasta")

# read in tree using ape
 tree = ape::read.tree("/home/david/ligases/raxml/alignment.fasta.raxml.supportFBP")
 rooted_tree = ape::root.phylo(tree, 
     outgroup = "D-ala--D-ala-ligase--Leuconostoc-mesenteroides")

# Plot with ggtree
p2 <- ggtree(rooted_tree) + 
    geom_tiplab(size=8, color="black", hjust = -0.05) +
    geom_tree(linetype = "solid", size = 2) +
    geom_nodelab(color = "grey40", size = 5, hjust=-0.3, size=5) +
    theme(legend.position="bottom",
        plot.margin = margin(1,12,1,1, "cm")) +
    coord_cartesian(clip = 'off') +
    hexpand(.5, direction = 1)
# save the tree to svg/png
ggsave(p2, file="/home/david/ligase_bootsrap2.svg", width = 12, height = 12)
ggsave(p2, file="/home/david/ligase_bootsrap2.png", width = 12, height = 12)




