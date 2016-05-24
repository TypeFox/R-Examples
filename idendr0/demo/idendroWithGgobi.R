## idendro + ggobi integration demo.
##

library(idendr0) # idendro

data(iris)

# perform hierarchical clustering analysis
hc <- hclust(dist(iris[, 1:4]))

# visualize clusters and heat map
idendro(hc, iris, hscale = 1, vscale = 1.2,
    ggobi = TRUE, ggobiGlyphType = 4, ggobiGlyphSize = 2)
