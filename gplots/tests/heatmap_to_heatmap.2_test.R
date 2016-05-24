library(gplots)
data(mtcars)

x <- as.matrix(mtcars)

## draws expected image
testHeatmap <- heatmap(x, Colv=NA, col=bluered(256), scale="column",
                       keep.dendro=TRUE)

# to prove this dendro is OK, redraw with same function:
heatmap(x, Colv=NA, col=bluered(256), scale="column",
        keep.dendro=TRUE, Rowv=testHeatmap$Rowv)

# but it doesn't work with heatmap.2()
heatmap.2(x, Colv=NA, col=bluered(256), scale="column",
          Rowv=testHeatmap$Rowv, dendrogram="row")
