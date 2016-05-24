require(phyloclim)

setwd("/Users/Stoffi/R_files/packages/phyloclim_archive/abellan")

phy <- read.nexus("grac_tree.NEX")
#phy <- fixTips(phy)
age <- branching.times(phy)

load("ovlap.RData")

# node 21 an 22 values of range overlap should be 0.0896 and 0.2913

# overlap <- sapply(names(age), nested.mean.overlap, phy = phy, 		olap = ovlap)

node <- 22
olap <- ovlap

plot(phy, no.margin = TRUE)
nodelabels(cex = 0.75)
tiplabels(cex = 0.75)
