### R code from vignette source 'gRbase-graphs.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: gRbase-graphs.Rnw:23-26
###################################################
require( gRbase )
prettyVersion <- packageDescription("gRbase")$Version
prettyDate <- format(Sys.Date())


###################################################
### code chunk number 2: gRbase-graphs.Rnw:70-74
###################################################
dir.create("fig")
oopt <- options()
options("digits"=4, "width"=80, "prompt"="R> ", "continue"="  ")
options(useFancyQuotes="UTF-8")


###################################################
### code chunk number 3: gRbase-graphs.Rnw:78-80
###################################################
library(gRbase)
library(graph)


###################################################
### code chunk number 4: gRbase-graphs.Rnw:197-200
###################################################
ug11 <- ug( ~a*b + b*c*d)
ug12 <- ug( c("a","b"), c("b","c","d") )
ug13 <- ug( list( c("a","b"), c("b","c","d") ) )


###################################################
### code chunk number 5: gRbase-graphs.Rnw:208-212
###################################################
ug11
nodes(ug11)
str( edges(ug11) )
plot( ug11 )


###################################################
### code chunk number 6: gRbase-graphs.Rnw:217-219
###################################################
ug11m <- ug( ~a*b + b*c*d, result="matrix")
ug11M <- ug( ~a*b + b*c*d, result="dgCMatrix")


###################################################
### code chunk number 7: gRbase-graphs.Rnw:230-233
###################################################
dag11 <- dag( ~a*b + b*c*d)
dag12 <- dag( c("a","b"), c("b","c","d") )
dag13 <- dag( list( c("a","b"), c("b","c","d") ) )


###################################################
### code chunk number 8: gRbase-graphs.Rnw:250-254
###################################################
dag11
nodes( dag11 )
str( edges( dag11 ) )
plot( dag11 )


###################################################
### code chunk number 9: gRbase-graphs.Rnw:260-262
###################################################
dag11m <- dag( ~a*b + b*c*d, result="matrix")
dag11M <- dag( ~a*b + b*c*d, result="dgCMatrix")


###################################################
### code chunk number 10: gRbase-graphs.Rnw:275-278
###################################################
d1.bi <- dag(~a:b + b:a)
edgemode( d1.bi )
str( edges(d1.bi) )


###################################################
### code chunk number 11: gRbase-graphs.Rnw:284-285
###################################################
d2.cyc <- dag(~a:b+b:c+c:a)


###################################################
### code chunk number 12: gRbase-graphs.Rnw:289-290
###################################################
par(mfrow=c(1,2)); plot(d1.bi); plot(d2.cyc)


###################################################
### code chunk number 13: gRbase-graphs.Rnw:296-297
###################################################
print( try( dag(~a:b+b:c+c:a, forceCheck=TRUE) ) )


###################################################
### code chunk number 14: gRbase-graphs.Rnw:388-391
###################################################
(mat <- as(ug11, "matrix"))
(Mat <- as(mat, "dgCMatrix"))
(NEL <- as(Mat, "graphNEL"))


###################################################
### code chunk number 15: gRbase-graphs.Rnw:396-399
###################################################
mat <- coerceGraph(ug11, "matrix")
Mat <- coerceGraph(ug11, "dgCMatrix")
NEL <- coerceGraph(mat, "graphNEL")


###################################################
### code chunk number 16: gRbase-graphs.Rnw:403-406
###################################################
mat <- graphNEL2M(ug11, result="matrix")
Mat <- graphNEL2M(ug11, result="dgCMatrix")
NEL <- M2graphNEL(mat)


###################################################
### code chunk number 17: gRbase-graphs.Rnw:411-417
###################################################
if( require(microbenchmark) ){
  microbenchmark(as(ug11, "matrix"), coerceGraph(ug11, "matrix"), graphNEL2M(ug11),
                 as(ug11, "dgCMatrix"), coerceGraph(ug11, "dgCMatrix"),
                 graphNEL2M(ug11, result="Matrix"),
                 as(mat, "graphNEL"), coerceGraph(mat, "graphNEL"), M2graphNEL(mat))
}


###################################################
### code chunk number 18: gRbase-graphs.Rnw:423-425
###################################################
str( edges(as(mat,"graphNEL")) )
str( M2adjList(mat) )


###################################################
### code chunk number 19: gRbase-graphs.Rnw:446-447
###################################################
dag11.mor <- moralize(dag11)


###################################################
### code chunk number 20: gRbase-graphs.Rnw:451-452
###################################################
par(mfrow=c(1,2)); plot(dag11); plot(dag11.mor)


###################################################
### code chunk number 21: gRbase-graphs.Rnw:459-461
###################################################
moralize( dag11m )
moralize( dag11, result="matrix" )


###################################################
### code chunk number 22: gRbase-graphs.Rnw:476-479
###################################################
topoSort(dag11)
topoSort(dag11m)
topoSort(dag11M)


###################################################
### code chunk number 23: gRbase-graphs.Rnw:486-487
###################################################
topoSort(dag(~a:b+b:c+c:a))


###################################################
### code chunk number 24: gRbase-graphs.Rnw:493-494
###################################################
topoSort( ug( ~a:b ) )


###################################################
### code chunk number 25: gRbase-graphs.Rnw:510-513
###################################################
str( getCliques(ug11) )
str( getCliques(ug11m) )
str( getCliques(ug11M) )


###################################################
### code chunk number 26: gRbase-graphs.Rnw:520-524
###################################################
if (require(microbenchmark)){
    microbenchmark(
        RBGL::maxClique( ug11 ), getCliques( ug11 ), getCliques( ug11m ),
        getCliques( ug11M ))}


###################################################
### code chunk number 27: gRbase-graphs.Rnw:541-544
###################################################
mcs(ug11)
mcs(ug11m)
mcs(ug11M)


###################################################
### code chunk number 28: gRbase-graphs.Rnw:550-551
###################################################
mcs(ug11, root=c("a","c"))


###################################################
### code chunk number 29: gRbase-graphs.Rnw:558-560
###################################################
mcs( dag11 )
mcs( as(dag11, "matrix") )


###################################################
### code chunk number 30: gRbase-graphs.Rnw:567-571
###################################################
ug11.nc <- ug( ~a:b:c + c:d + d:e + a:e + f:g )
mcs(ug11.nc)
(tug11.nc  <- triangulate(ug11.nc))
mcs(tug11.nc)


###################################################
### code chunk number 31: gRbase-graphs.Rnw:575-576
###################################################
par(mfrow=c(1,2)); plot(ug11.nc); plot(tug11.nc)


###################################################
### code chunk number 32: gRbase-graphs.Rnw:594-595
###################################################
rp <- rip(tug11.nc); rp


###################################################
### code chunk number 33: gRbase-graphs.Rnw:599-600
###################################################
plot( rp )


###################################################
### code chunk number 34: gRbase-graphs.Rnw:623-625
###################################################
g1 <- ug(~a:b+b:c+c:d+d:e+e:f+a:f+b:e)
g1mt <- minimalTriang(g1) # A minimal triangulation


###################################################
### code chunk number 35: gRbase-graphs.Rnw:629-630
###################################################
par(mfrow = c(1,2)); plot(g1); plot(g1mt)


###################################################
### code chunk number 36: gRbase-graphs.Rnw:636-638
###################################################
g2 <- ug(~a:b:e:f+b:c:d:e)
g1mt2 <- minimalTriang(g1, tobject=g2)


###################################################
### code chunk number 37: gRbase-graphs.Rnw:642-643
###################################################
par(mfrow = c(1,2)); plot(g2); plot(g1mt2)


###################################################
### code chunk number 38: gRbase-graphs.Rnw:648-649
###################################################
mm <- mpd( g1 ); mm


###################################################
### code chunk number 39: gRbase-graphs.Rnw:653-656
###################################################
par(mfrow = c(1,2))
plot(subGraph(mm$cliques[[1]], g1))
plot(subGraph(mm$cliques[[2]], g1))


###################################################
### code chunk number 40: gRbase-graphs.Rnw:674-680
###################################################
if(require(microbenchmark)){
    microbenchmark(
        RBGL::maxClique(ug11),
        getCliques(ug11),
        getCliques(ug11m),
        getCliques(ug11M)  ) }


###################################################
### code chunk number 41: gRbase-graphs.Rnw:696-712
###################################################
V <- 1:300
M <- 1:10
## Sparse graph
##
g1 <- randomGraph(V, M, 0.05)
length(edgeList(g1))
s<-c(NEL=object.size(g1), dense=object.size(as(g1, "matrix")),
     sparse=object.size(as(g1, "dgCMatrix")))
s/max(s)
## More dense graph
##
g1 <- randomGraph(V, M, 0.5)
length(edgeList(g1))
s <- c(NEL=object.size(g1), dense=object.size(as(g1, "matrix")),
       sparse=object.size(as(g1, "dgCMatrix")))
s/max(s)


###################################################
### code chunk number 42: gRbase-graphs.Rnw:730-731
###################################################
args(querygraph)


###################################################
### code chunk number 43: gRbase-graphs.Rnw:738-740
###################################################
#rm(print.list)
options("width"=85)


