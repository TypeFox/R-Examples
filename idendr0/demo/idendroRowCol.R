## Demo on clustering both rows and columns of a data set.
##

library(idendr0) # idendro

data(iris)

# Find similar variables (i.e. features measured on the Iris flowers:
# sepal/petal widths/lengths) ...
hx <- hclust(dist(t(scale(iris[, 1:4]))))
# ... and display the similarities in a dendrogram.
idendro(hx, t(iris[, 1:4]), heatmapRelSize=.7, doScaleHeatmapByRows=TRUE, mai=par('mai')*c(1,1,1,2))
# (Petal width and length being most similar.)
