library(phylobase)
library(ape)

data(geospiza)
g1 <- as(geospiza,"phylo4")
g2 <- geospiza

par(mfrow=c(1,2))
plot(g1, show.node.label=TRUE)
## be careful with this: works if par("fin")=c(5.56,6.77)
##                       fails if par("fin")=c(4.87,6.77)
##try(plot(g2,show.node.label=TRUE),silent=TRUE)
## Here, R was complaining about a lack of room to plot data
## so nothing abnormal. -- TJ
plot(g2, show.node.label=TRUE)


## commented out since phylog objects are deprecated anyway
## g2B <- as(extractTree(g2), "phylog")
##  Note the numbering differences!

## round trip
g2C <- as(read.tree(text=write.tree(as(g1, "phylo"))), "phylo4")
## comes back in same order
try(plot(g1, show.node.label=TRUE))
try(plot(g2C, show.node.label=TRUE))

g3 = subset(g2, tips.exclude=c("fuliginosa", "fortis", "magnirostris",
                 "conirostris", "scandens"))
plot(extractTree(g3))  ## phylo4
plot(g3)


## Playing with new ways of plotting

if (FALSE) {
if(require(MASS)){
    dist1 <- cophenetic.phylo(as(g2, "phylo"))
    mdspos <- isoMDS(dist1)$points
    par(mfrow=c(2, 2))
    plot(g1)
    ## plot(mdspos,type="n")
    ## text(mdspos[,1],mdspos[,2],abbreviate(rownames(mdspos)))
    ## cmdpos <- cmdscale(dist1)
    ## plot(cmdpos,type="n")
    ## text(cmdpos[,1],cmdpos[,2],abbreviate(rownames(mdspos)))
}
## never mind, I don't know how to construct a useful
##  2D color space anyway ...
}

treePlot(g2,plot.at.tip=TRUE,tip.plot.fun=
         function(x,...) {
           grid::grid.points(seq(along=x),x)})
