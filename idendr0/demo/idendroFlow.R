## Visualization of flow cytometry data in dendrogram and heat
## map, with two scatter plots integrated.
##

library(flowStats) # ITN
library(RColorBrewer) # brewer.pal
library(idendr0) # idendro

data(ITN)

# get data matrix
x <- exprs(ITN$sample03[, 1:7])
# preprocess the data
x[, 3:7] <- log10(x[, 3:7])
x <- scale(x)

# perform HCA
cat('computing HCA (this takes a while)\n')
hx <- hclust(dist(x), method = 'average')

# setup scatter plots of selected data projections
dev.new()
par(ask = FALSE)
sp1 <- dev.cur()
dev.new()
par(ask = FALSE)
sp2 <- dev.cur()

clusterColors <- brewer.pal(12, "Paired")

colorizeCallback <- function(clr) {
    d<-dev.cur()
    dev.set(sp1)
    plot(x[,"CD3"], x[,"HLADr"], pch = 19, col = c('black',clusterColors)[clr + 1])
    dev.set(sp2)
    plot(x[,"CD8"], x[,"CD4"], pch = 19, col = c('black',clusterColors)[clr + 1])
    dev.set(d)
}

# produce the scatter plots by invoking the callback function with no specific colors
colorizeCallback(0)

# plot dendrogram + heat map
cat('plotting dendrogram\n')
idendro(hx, x,
    heatmapColors = colorRampPalette(c("purple4", "blue3", "blue3", "grey",
        "grey", "orangered", "orangered", "red"))(15),
    clusterColors = clusterColors,
    colorizeCallback = colorizeCallback)

