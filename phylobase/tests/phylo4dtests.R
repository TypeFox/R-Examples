library(phylobase)
library(ape)
tree.phylo <- read.tree(text="(((A,B)C,D),E);")  #only one node is labelled
tree <- as(tree.phylo, "phylo4")

tree.phylo2 <- read.tree(text="(((A,B)C,D)F,E)G;")  # all nodes labelled
tree2 <- as(tree.phylo2, "phylo4")

tip.data <- data.frame(size=c(1, 2, 3, 4))
rownames(tip.data) <- c("A", "B", "E", "D")

treed <- phylo4d(tree, tip.data)
dat2 <- data.frame(size=c(0,1,2), row.names=c("G", "F", "C"))

try(phylo4d(tree, node.data=dat2), silent = TRUE)  # error, cannot match data because no node labels on tree
phylo4d(tree2, node.data=dat2) -> treed2  # OK tree labelled; has node data, no tip data 

plot(treed2) # works with a warning about no tip data to plot
tipData(treed2, empty.columns=FALSE) #returns empty 4-row data.frame

phylo4d(tree2, tip.data=tip.data, node.data=dat2) -> treed3 #node+tip data

plot(treed3)  # works
tipData(treed3)  #works, but returns tips only
tdata(treed3, "all")

print(tree)
print(treed)

