require(phyloclim)

# setwd("/Users/Stoffi/R_files/packages/phyloclim_archive/abellan")
# phy <- read.nexus("grac_tree.NEX")

## load PNOs for Oxalis sect. Palmatifoliae ...
data(PNO)
## ... and calculate niche overlap for annual mean temperature
no <- niche.overlap(PNO$AnnualMeanTemperature)

## load chronogram for section Palmatifoliae
data(tree); plot(tree)

x <- age.range.correlation(phy = tree, overlap = no, n = 100)

phy <- tree; n = 10; tri = "upper"; overlap <- no