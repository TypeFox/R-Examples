require("GMD")          # load library
data(mtcars)            # load data
x  <- as.matrix(mtcars) # data as a matrix

dev.new(width=10,height=8)
heatmap.3(x)                               # default, with reordering and dendrogram
heatmap.3(x, Rowv=FALSE, Colv=FALSE)       # no reordering and no dendrogram
heatmap.3(x, dendrogram="none")            # reordering without dendrogram
heatmap.3(x, dendrogram="row")             # row dendrogram with row (and col) reordering
heatmap.3(x, dendrogram="row", Colv=FALSE) # row dendrogram with only row reordering
heatmap.3(x, dendrogram="col")             # col dendrogram
heatmap.3(x, dendrogram="col", Rowv=FALSE) # col dendrogram with only col reordering
heatmapOut <-
  heatmap.3(x, scale="column")             # sacled by column
names(heatmapOut)                          # view the list that is returned
heatmap.3(x, scale="column", x.center=0)   # colors centered around 0
heatmap.3(x, scale="column",trace="column")# trun "trace" on

## coloring cars (row observations) by brand
brands <- sapply(rownames(x), function(e) strsplit(e," ")[[1]][1]) 
names(brands) <- c()
brands.index <- as.numeric(as.factor(brands))
RowIndividualColors <- rainbow(max(brands.index))[brands.index]
heatmap.3(x, scale="column", RowIndividualColors=RowIndividualColors)

## coloring attributes (column features) randomly (just for a test :)
heatmap.3(x, scale="column", ColIndividualColors=rainbow(ncol(x)))

## add a single plot for all row individuals 
dev.new(width=12,height=8)
expr1 <- list(quote(plot(row.data[rowInd,"hp"],1:nrow(row.data),xlab="hp",ylab="",
                         main="Gross horsepower",yaxt="n")),
              quote(axis(2,1:nrow(row.data),rownames(row.data)[rowInd],las=2)))
heatmap.3(x, scale="column", plot.row.individuals=TRUE, row.data=x,
          plot.row.individuals.list=list(expr1))

## add a single plot for all col individuals 
dev.new(width=12,height=8)
expr2 <- list(quote(plot(colMeans(col.data)[colInd], xlab="", ylab="Mean",xaxt="n",
                         main="Mean features",cex=1,pch=19)),
              quote(axis(1,1:ncol(col.data),colnames(col.data)[colInd],las=2)))
heatmap.3(x, scale="column", plot.col.individuals=TRUE, col.data=x,
          plot.col.individuals.list=list(expr2))

## add another single plot for all col individuals
dev.new(width=12,height=8)
expr3 <- list(quote(op <- par(mar = par("mar")*1.5)),
              quote(mytmp.data <- apply(col.data,2,function(e) e/sum(e))),
              quote(boxplot(mytmp.data[,colInd], xlab="", ylab="Value",
                         main="Boxplot of normalized column features")),
              quote(par(op)))
heatmap.3(x, scale="column", plot.col.individuals=TRUE, col.data=x,
          plot.col.individuals.list=list(expr3))
