## idendro + scatter plot + parallel coordinate plot integration demo.
##

library(idendr0) # idendro

data(iris)

# perform hierarchical clustering analysis
hc <- hclust(dist(iris[, 1:4]))

# setup a scatter plot
dev.new()
par(ask = FALSE)
scatterDevId <- dev.cur()

# setup a parallel coordinate plot
if (require(MASS)) {
    dev.new()
    par(ask = FALSE)
    parcoordDevId <- dev.cur()
    iris.numeric <- iris
    iris.numeric$Species <- as.numeric(iris.numeric$Species)
} else {
    parcoordDevId <- NA
}

colorizeCallback <- function(clr) {
    clusterColors <- c('black', 'red', 'green', 'blue', 'yellow', 'magenta',
        'cyan', 'darkred', 'darkgreen', 'purple', 'darkcyan')

    # color the scatter plot
    dev.set(scatterDevId)
    opar <- par(ask = FALSE, mai = par('mai') * c(1, 1, .5, 1))
    with(iris, plot(Sepal.Length, Sepal.Width, col = clusterColors[clr + 1], ty = 'p', pch=19))
    par(opar)

    # color the parallel coordinate plot
    if (!is.na(parcoordDevId)) {
        dev.set(parcoordDevId)
        opar <- par(ask = FALSE, mai = par('mai') * c(.5, 1, 1, 1))
        parcoord(iris.numeric, col = clusterColors[clr + 1])
        par(opar)
    }
}

# produce the scatter plot and the parallel coordinate plot
# by invoking the callback function with no specific colors
colorizeCallback(0)
# visualize clusters and heat map
idendro(hc, iris, colorizeCallback = colorizeCallback, hscale = 1.3, vscale = 1.3)
