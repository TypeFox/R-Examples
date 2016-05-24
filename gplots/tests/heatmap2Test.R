library(gplots)

data(mtcars)
x  <- as.matrix(mtcars)
rc <- rainbow(nrow(x), start=0, end=.3)
cc <- rainbow(ncol(x), start=0, end=.3)

heatmap.2(x)  ## default - dendrogram plotted and reordering done.
heatmap.2(x, dendrogram="none") ##  no dendrogram plotted, but reordering done.
heatmap.2(x, dendrogram="row") ## row dendrogram plotted and row reordering done.
heatmap.2(x, dendrogram="col") ## col dendrogram plotted and col reordering done.

heatmap.2(x, keysize=2)  ## default - dendrogram plotted and reordering done.


heatmap.2(x, Colv=FALSE)
heatmap.2(x, Colv=NULL)

heatmap.2(x, Rowv=FALSE)
heatmap.2(x, Rowv=NULL)

heatmap.2(x, Rowv=FALSE, dendrogram="both") ## generate warning!
heatmap.2(x, Rowv=NULL, dendrogram="both")  ## generate warning!
heatmap.2(x, Colv=FALSE, dendrogram="both") ## generate warning!



heatmap.2(x,dendrogram="both", Rowv=FALSE) ## row and col reordering done, but only col dendrogram plotted 
heatmap.2(x,dendrogram="row",  Rowv=FALSE) 


heatmap.2(x,dendrogram="col",  Rowv=FALSE)
heatmap.2(x,dendrogram="none", Rowv=FALSE) 


heatmap.2(x,dendrogram="both", Colv=FALSE) ## row and col reordering done, but only row dendrogram plotted
heatmap.2(x,dendrogram="row",  Colv=FALSE)
heatmap.2(x,dendrogram="col",  Colv=FALSE)
heatmap.2(x,dendrogram="none", Colv=FALSE) 

## above statement is same as heatmap.2(x,Colv=FALSE)
heatmap.2(x,dendrogram="col", Colv=FALSE) ## col reordering done, but no dendrogram plotted 


hv <- heatmap.2(x, col=cm.colors(256), scale="column", 
	       RowSideColors=rc, ColSideColors=cc, margin=c(5, 10), 
	       xlab="specification variables", ylab= "Car Models", 
	       main="heatmap(<Mtcars data>, \n        ..., scale=\"column\")", 
               tracecol="green", density="density")

str(hv) # the two re-ordering index vectors

## Test row and column separators 
heatmap.2(x,dendrogram="col", colsep=1) ## one colsep
heatmap.2(x,dendrogram="col", rowsep=1) ## one rowsep

heatmap.2(x,dendrogram="col", colsep=1:2) ## two colseps
heatmap.2(x,dendrogram="col", rowsep=1:2) ## two rowseps
