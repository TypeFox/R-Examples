## Visualization of spectroscopy data in dendrogram and heat map with
## a scatter plot and parallel coordinate plot integrated.
##

library(hyperSpec) # chondro data
library(MASS) # parcoord
library(idendr0) # idendro

# preprocess chondro data
cat('preprocessing the chondro data (it takes a while)\n')
chondro <- spc.loess(chondro, newx = seq (602, 1800, by = 4))
chondro <- chondro - spc.fit.poly.below(chondro)
chondro <- sweep(chondro, 1, rowMeans(chondro), "/")
overall.composition <- quantile(chondro, 0.05)
chondro <- sweep(chondro, 2, overall.composition, "-")

# HCA
dst <- dist(chondro)
dndr <- hclust(dst, method = "ward")

# data frame construction
df.chondro <- as.wide.df(chondro)
colnames(df.chondro)[-(1:3)] <- paste("wl", colnames(df.chondro)[-(1:3)], sep = ".")
names <- as.character(wl(chondro))
names[wl(chondro) %% 50 != 0] <- ""

# setup a scatter plot
dev.new()
par(ask = FALSE)
p1 <- dev.cur()

# setup a parallel coordinate plot
dev.new()
par(ask = FALSE)
p2 <- dev.cur()

colorizeCallback <- function(clr) {
    clusterColors <- c('black', 'red', 'green', 'blue', 'yellow', 'magenta',
        'cyan', 'darkred', 'darkgreen', 'purple', 'darkcyan')

    # color the scatter plot
    dev.set(p1)
    with(df.chondro, plot(x, y, pch = 19, col = clusterColors[clr + 1]))

    dev.set(p2)
    parcoord(df.chondro[, !colnames(df.chondro) %in% 'clusters'], col = clusterColors[clr + 1])
}

# produce the two plots by invoking the callback function with no specific colors
colorizeCallback(0)

# dendrogram
idendro(dndr, df.chondro, heatmapRelSize = 0.75, heatmapColors = alois.palette(25),
    colorizeCallback = colorizeCallback)
