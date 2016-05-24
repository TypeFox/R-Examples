require("GMD")      # load library
data(chipseq_mES)   # load data
data(chipseq_hCD4T) # load data

## pairwise distance and alignment based on GMD metric
plot(gmdm(chipseq_mES,sliding=FALSE))

## clustering on spatial distributions of histone modifications
x <- gmdm(chipseq_hCD4T,sliding=FALSE,resolution=10)
heatmap.3(x,revC=TRUE) 

## Determine the number of clusters by "Elbow" criterion
main <- "Heatmap of ChIP-seq data (human CD4+ T cells)"
dist.obj <- gmdm2dist(x)
css.multi.obj <- css.hclust(dist.obj,hclust(dist.obj))
elbow.obj <- elbow.batch(css.multi.obj,ev.thres=0.90,inc.thres=0.05)
heatmap.3(dist.obj, main=main, revC=TRUE, kr=elbow.obj$k, kc=elbow.obj$k)

## more strict threshold 
elbow.obj <- elbow.batch(css.multi.obj,ev.thres=0.75,inc.thres=0.1)
heatmap.3(dist.obj, main=main, revC=TRUE, kr=elbow.obj$k, kc=elbow.obj$k)

## side plots 
normalizeVector <- function(v){v/sum(v)} # a function to normalize a vector
dev.new(width=12,height=8)

## summary of row clusters
expr1 <- list(quote(op <- par(mar = par("mar")*2)),
              quote(plot(mhist.summary(as.mhist(i.x)),if.plot.new=FALSE)),
              quote(par(op))
              )

## summary of row clustering
expr2 <- list(quote(tmp.clusters <- cutree(hclust(dist.row),k=kr)),
              quote(tmp.css <- css(dist.row,tmp.clusters)),
              quote(print(tmp.css)),
              quote(tmp.wev <- tmp.css$wss/tmp.css$tss),
              quote(names(tmp.wev) <- as.character(unique(tmp.clusters))),
              quote(tmp.wev <- tmp.wev[order(unique(tmp.clusters))]),
              quote(barplot(tmp.wev,main="Cluster Explained Variance", xlab="Cluster",
                            ylab="EV",col="white",border="black",
                            ylim=c(0,max(tmp.wev)*1.1),cex.main=1)))
expr3 <- list(quote(op <- par(mar = par("mar")*2)),
              quote(plot.elbow(css.multi.obj,elbow.obj,if.plot.new=FALSE)),
              quote(par(op))
              )

heatmap.3(dist.obj, main=main, cex.main=1.25, revC=TRUE, kr=elbow.obj$k, kc=elbow.obj$k,
          keysize=1,mapsize=4.5,row.data=lapply(chipseq_hCD4T,normalizeVector),
          plot.row.clusters=TRUE,plot.row.clusters.list=list(expr1),
          plot.row.clustering=TRUE,plot.row.clustering.list=list(expr2,expr3))

