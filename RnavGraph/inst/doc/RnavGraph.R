### R code from vignette source 'RnavGraph.Rnw'

###################################################
### code chunk number 1: RnavGraph.Rnw:58-60 (eval = FALSE)
###################################################
## library(tcltk)
## .Tcl('set tcl_version')


###################################################
### code chunk number 2: RnavGraph.Rnw:65-69 (eval = FALSE)
###################################################
## source("http://www.bioconductor.org/biocLite.R")
## biocLite("graph"); biocLite("RBGL"); biocLite("Rgraphviz"); biocLite("RDRToolbox")
## install.packages(c("PairViz", "scagnostics", "rgl", "MASS",
##       "hexbin", "vegan", "png"), dependencies = TRUE)


###################################################
### code chunk number 3: RnavGraph.Rnw:73-75 (eval = FALSE)
###################################################
## install.packages("RnavGraph")
## install.packages("RnavGraphImageData")


###################################################
### code chunk number 4: RnavGraph.Rnw:88-89
###################################################
library(RnavGraph)


###################################################
### code chunk number 5: RnavGraph.Rnw:92-93
###################################################
ls("package:RnavGraph")


###################################################
### code chunk number 6: RnavGraph.Rnw:98-99 (eval = FALSE)
###################################################
## demo(package = "RnavGraph")


###################################################
### code chunk number 7: RnavGraph.Rnw:102-103 (eval = FALSE)
###################################################
## system.file("demo", package="RnavGraph")


###################################################
### code chunk number 8: RnavGraph.Rnw:112-117
###################################################
library(RnavGraph)
ng.iris <- ng_data(name = "iris", data = iris[,1:4],
                   shortnames = c('s.L', 's.W', 'p.L', 'p.W'),
                   group = iris$Species,
                   labels = substr(iris$Species,1,2))


###################################################
### code chunk number 9: RnavGraph.Rnw:123-124 (eval = FALSE)
###################################################
## navGraph(ng.iris)


###################################################
### code chunk number 10: RnavGraph.Rnw:345-346 (eval = FALSE)
###################################################
## demo("ng_2d_iris", package = "RnavGraph")


###################################################
### code chunk number 11: RnavGraph.Rnw:351-355
###################################################
ng.iris <- ng_data(name = "iris", data = iris[,1:4],
                   shortnames = c('s.L', 's.W', 'p.L', 'p.W'),
                   group = iris$Species,
                   labels = substr(iris$Species,1,2))


###################################################
### code chunk number 12: RnavGraph.Rnw:359-360
###################################################
V <- shortnames(ng.iris)


###################################################
### code chunk number 13: RnavGraph.Rnw:363-366
###################################################
G <- completegraph(V)
LG <- linegraph(G)
LGnot <- complement(LG)


###################################################
### code chunk number 14: RnavGraph.Rnw:369-370
###################################################
class(G)


###################################################
### code chunk number 15: RnavGraph.Rnw:373-375
###################################################
ng.lg <- ng_graph(name = '3D Transition', graph = LG, layout = 'circle')
ng.lgnot <- ng_graph(name = '4D Transition', graph = LGnot, layout = 'circle')


###################################################
### code chunk number 16: RnavGraph.Rnw:380-381
###################################################
viz3dTransition <- ng_2d(ng.iris,ng.lg)


###################################################
### code chunk number 17: RnavGraph.Rnw:384-385
###################################################
viz4dTransition  <- ng_2d(ng.iris,ng.lgnot)


###################################################
### code chunk number 18: RnavGraph.Rnw:388-390
###################################################
viz <- list(viz3dTransition, viz4dTransition)	
graphs <- list(ng.lg, ng.lgnot)


###################################################
### code chunk number 19: RnavGraph.Rnw:393-394 (eval = FALSE)
###################################################
## nav <- navGraph(graph = graphs, data = ng.iris, viz = viz)


###################################################
### code chunk number 20: RnavGraph.Rnw:411-412
###################################################
ng.iris <- ng_data(name = "iris", data = iris[,1:4])


###################################################
### code chunk number 21: RnavGraph.Rnw:415-416
###################################################
ng.iris


###################################################
### code chunk number 22: RnavGraph.Rnw:421-423
###################################################
names(ng.iris)  ## get variable names
names(ng.iris) <- c("SepL", "SepL",  "PetL", "PetW")


###################################################
### code chunk number 23: RnavGraph.Rnw:428-430
###################################################
shortnames(ng.iris)
shortnames(ng.iris) <- c('s.L', 's.W', 'p.L', 'p.W')


###################################################
### code chunk number 24: RnavGraph.Rnw:435-436
###################################################
print(ng.iris)


###################################################
### code chunk number 25: RnavGraph.Rnw:441-446
###################################################
ng_get(ng.iris) ## See what can be accessed
ng_get(ng.iris,"data")[1:3,]
ng_get(ng.iris,"name")
ng_get(ng.iris,"group")
ng_get(ng.iris,"labels")


###################################################
### code chunk number 26: RnavGraph.Rnw:451-453
###################################################
ng_set(ng.iris)
ng_set(ng.iris,"labels") <- as.character(iris$Species)


###################################################
### code chunk number 27: RnavGraph.Rnw:460-464
###################################################
x <- rep(1:30, each = 30)
y <- rep(1:30, 30)
ng.test <- ng_data(name = "test", data = data.frame(x = x, y = y, z = 1:900),
		group = 1:900)


###################################################
### code chunk number 28: RnavGraph.Rnw:469-470 (eval = FALSE)
###################################################
## nav <- navGraph(ng.test)


###################################################
### code chunk number 29: RnavGraph.Rnw:484-489 (eval = FALSE)
###################################################
## x <- rep(1:3, each = 3)
## y <- rep(1:3, 3)
## ng.test <- ng_data(name = "test2", data = data.frame(x = x, y = y, z = 1:9),
## 		group = c('red','red','green','blue','blue','blue','blue','orange','orange'))
## nav <- navGraph(ng.test)


###################################################
### code chunk number 30: RnavGraph.Rnw:508-511 (eval = FALSE)
###################################################
## vignette(package="graph")
## vignette(package="RBGL")
## vignette(package="Rgraphviz")


###################################################
### code chunk number 31: RnavGraph.Rnw:516-521
###################################################
adjM <- matrix(c(0,4,1,0,2,0,3,2,2,2,0,0,0,2,0,0), ncol = 4)
rownames(adjM) <- c('A','B','C','D')
colnames(adjM) <- c('A','B','C','D')
G <- new("graphAM", adjMat = adjM, edgemode = "directed")
G


###################################################
### code chunk number 32: RnavGraph.Rnw:526-528
###################################################
library(Rgraphviz)
plot(G)


###################################################
### code chunk number 33: RnavGraph.Rnw:535-538
###################################################
V <-  c('s.L', 's.W', 'p.L', 'p.W')
G <- completegraph(V)
plot(G, "neato")


###################################################
### code chunk number 34: RnavGraph.Rnw:545-551
###################################################
from <- c("A","A","C","C")
to   <- c("B","C","B","D")
ftEmat <- cbind(from,to)

G <- newgraph(nodeNames = LETTERS[1:5], mat = ftEmat)
plot(G, "neato")


###################################################
### code chunk number 35: RnavGraph.Rnw:557-559
###################################################
G <- newgraph(nodeNames = LETTERS[1:5], mat = ftEmat, directed = TRUE)
plot(G, "neato")


###################################################
### code chunk number 36: RnavGraph.Rnw:565-570
###################################################
from <- c(1,1,3,3)
to   <- c(2,3,2,4)
ftEmat <- cbind(from,to)
G <- newgraph(nodeNames = LETTERS[1:5], mat = ftEmat)
plot(G, "neato")


###################################################
### code chunk number 37: RnavGraph.Rnw:575-578
###################################################
weights <- c(2,1,3,4)
G <- newgraph(nodeNames = LETTERS[1:5], mat = ftEmat, weights = weights)
edgeData(G, attr = "weight")


###################################################
### code chunk number 38: RnavGraph.Rnw:582-586
###################################################
G <- newgraph(nodeNames = LETTERS[1:5], mat = ftEmat, weights = weights, directed = TRUE)
edgeData(G, attr = "weight")
edgeData(G, from = "A", to = "B", attr = "weight")
edgeData(G, from = "A", to = "B", attr = "weight") <- 8


###################################################
### code chunk number 39: RnavGraph.Rnw:591-595
###################################################
adjM <- matrix(c(0,1,1,0,1,0,1,1,1,1,0,0,0,1,0,0), ncol = 4)
all(adjM == t(adjM)) ## is symmetric (undirected)
G <- newgraph(nodeNames = V, mat= adjM, isAdjacency=TRUE)
plot(G,"neato")


###################################################
### code chunk number 40: RnavGraph.Rnw:599-601
###################################################
nodes(G)
edges(G)


###################################################
### code chunk number 41: RnavGraph.Rnw:604-608
###################################################
adjM <-     matrix(c(0,0,1,0,1,0,1,1,0,0,0,0,0,1,0,0), ncol = 4)
weightsM <- matrix(c(0,0,5,0,2,0,1,3,0,0,0,0,0,7,0,0), ncol = 4)
G <- newgraph(nodeNames = V, mat= adjM, weights = weightsM, directed = TRUE, isAdjacency=TRUE)
edgeData(G, attr = "weight")


###################################################
### code chunk number 42: RnavGraph.Rnw:611-612
###################################################
ftM2adjM(ftEmat)


###################################################
### code chunk number 43: RnavGraph.Rnw:618-621
###################################################
G <- completegraph(V)
LG <- linegraph(G, sep = '::')
nodes(LG)


###################################################
### code chunk number 44: RnavGraph.Rnw:633-635
###################################################
ng.LG <- ng_graph(name = "3D Transition", graph = LG, sep = '++', layout = "circle")
ng.LG  


###################################################
### code chunk number 45: RnavGraph.Rnw:638-641
###################################################
LGnot <- complement(LG)
ng.LGnot <- ng_graph(name = "4D Transition", 
                     graph = LGnot, sep = "::", layout = "circle")


###################################################
### code chunk number 46: RnavGraph.Rnw:646-649
###################################################
par(mfrow = c(1,2))
plot(ng.LG)
plot(ng.LGnot)


###################################################
### code chunk number 47: RnavGraph.Rnw:654-656
###################################################
ng_get(ng.LG)
ng_set(ng.LG, "name") <- "3d transition graph"


###################################################
### code chunk number 48: RnavGraph.Rnw:761-763
###################################################
viz1 <- ng_2d(data = ng.iris, graph = ng.lg)
viz1


###################################################
### code chunk number 49: RnavGraph.Rnw:766-767
###################################################
viz2 <- ng_2d(data = ng.iris, graph = ng.lgnot)


###################################################
### code chunk number 50: RnavGraph.Rnw:772-773
###################################################
library(RnavGraphImageData)


###################################################
### code chunk number 51: RnavGraph.Rnw:776-777 (eval = FALSE)
###################################################
## data(package='RnavGraphImageData')


###################################################
### code chunk number 52: RnavGraph.Rnw:782-784
###################################################
data(digits)
dim(digits)


###################################################
### code chunk number 53: RnavGraph.Rnw:787-788 (eval = FALSE)
###################################################
## help("digits")


###################################################
### code chunk number 54: RnavGraph.Rnw:791-792 (eval = FALSE)
###################################################
## matrix(digits[,7*1100+1],ncol = 16, byrow=FALSE)


###################################################
### code chunk number 55: RnavGraph.Rnw:795-799
###################################################
w <- options("width")$width
options(width=120)
matrix(digits[,7*1100+1],ncol = 16, byrow=FALSE)
options(width=w)


###################################################
### code chunk number 56: RnavGraph.Rnw:803-805
###################################################
sel <- sample(x=1:11000,size = 600)
p.digits <- digits[,sel]


###################################################
### code chunk number 57: RnavGraph.Rnw:809-812
###################################################
ng.i.digits <- ng_image_array_gray('USPS Handwritten Digits',
                                   p.digits,16,16,invert = TRUE,
                                   img_in_row = FALSE)


###################################################
### code chunk number 58: RnavGraph.Rnw:816-817 (eval = FALSE)
###################################################
## ng.i.digits


###################################################
### code chunk number 59: RnavGraph.Rnw:824-825
###################################################
library(vegan)


###################################################
### code chunk number 60: RnavGraph.Rnw:828-829 (eval = FALSE)
###################################################
## p.digitsT <- t(p.digits)


###################################################
### code chunk number 61: RnavGraph.Rnw:832-834 (eval = FALSE)
###################################################
## dise <- vegdist(p.digitsT, method="euclidean")
## ord <- isomap(dise,k = 8, ndim=6, fragmentedOK = TRUE)


###################################################
### code chunk number 62: RnavGraph.Rnw:837-844 (eval = FALSE)
###################################################
## digits_group <- rep(c(1:9,0), each = 1100)
## 
## ng.iso.digits <- ng_data(name = "ISO_digits",
## 		data = data.frame(ord$points),
## 		shortnames = paste('i',1:6, sep = ''),
## 		group = digits_group[sel],
## 		labels = as.character(digits_group[sel]))


###################################################
### code chunk number 63: RnavGraph.Rnw:847-853 (eval = FALSE)
###################################################
## V <- shortnames(ng.iso.digits)
## G <- completegraph(V)
## LG <- linegraph(G)
## LGnot <- complement(LG)
## ng.LG <- ng_graph(name = "3D Transition", graph = LG)
## ng.LGnot <- ng_graph(name = "4D Transition", graph = LGnot)


###################################################
### code chunk number 64: RnavGraph.Rnw:857-859 (eval = FALSE)
###################################################
## vizDigits1 <- ng_2d(data = ng.iso.digits, graph = ng.LG, images = ng.i.digits)
## vizDigits2 <- ng_2d(data = ng.iso.digits, graph = ng.LGnot, images = ng.i.digits)


###################################################
### code chunk number 65: RnavGraph.Rnw:862-863 (eval = FALSE)
###################################################
## nav <- navGraph(data = ng.iso.digits, graph = list(ng.LG, ng.LGnot), viz = list(vizDigits1, vizDigits2))


###################################################
### code chunk number 66: RnavGraph.Rnw:877-878 (eval = FALSE)
###################################################
## shell("convert image.png image.jpg")


###################################################
### code chunk number 67: RnavGraph.Rnw:906-908 (eval = FALSE)
###################################################
## imgPath <- "path to folder with png files"
## aloi_images <- list.files(path=imgPath, full.names=TRUE)


###################################################
### code chunk number 68: RnavGraph.Rnw:913-915 (eval = FALSE)
###################################################
## sel <- sample(1:length(aloi_images),replace=FALSE)
## p.aloi_images <- aloi_images[sel]


###################################################
### code chunk number 69: RnavGraph.Rnw:921-922 (eval = FALSE)
###################################################
## ng.i.objects <- ng_image_files(name="ALOI_Objects_Images", path=p.aloi_images)


###################################################
### code chunk number 70: RnavGraph.Rnw:925-926 (eval = FALSE)
###################################################
## ng.i.objects


###################################################
### code chunk number 71: RnavGraph.Rnw:930-938 (eval = FALSE)
###################################################
## library(png)
## imgData <- t(sapply(p.aloi_images, FUN=function(path){
##   x <- readPNG(path)
##   r <- sum(x[,,1])
##   g <- sum(x[,,2])
##   b <- sum(x[,,3])
##   return(c(r,g,b))
## }))


###################################################
### code chunk number 72: RnavGraph.Rnw:941-960 (eval = FALSE)
###################################################
## ng.iso.objects <- ng_data(name = "Aloi_Objects",
## 		data = data.frame(imgData),
## 		shortnames = paste('i',1:3, sep = ''))
## 
## ## 3d and 4d transition graph Graphs
## V <- shortnames(ng.iso.objects)
## G <- completegraph(V)
## LG <- linegraph(G)
## LGnot <- complement(LG)
## ng.LG <- ng_graph(name = "3D Transition", graph = LG)
## ng.LGnot <- ng_graph(name = "4D Transition", graph = LGnot)
## 
## ## visualization instructions
## vizObjects1 <- ng_2d(data = ng.iso.objects, graph = ng.LG, images = ng.i.objects)
## vizObjects2 <- ng_2d(data = ng.iso.objects, graph = ng.LGnot, images = ng.i.objects)
## 
## ## start navGraph 
## nav <- navGraph(data = ng.iso.objects, graph = list(ng.LG, ng.LGnot),
##                 viz = list(vizObjects1, vizObjects2))


###################################################
### code chunk number 73: RnavGraph.Rnw:1062-1067
###################################################
vizGlyph1 <- ng_2d(data = ng.iris, graph=ng.lg,
                   glyph=names(ng.iris)[c(1,2,3,4,3,2,3,4,1)])
vizGlyph2 <- ng_2d(data = ng.iris, graph=ng.lgnot, 
                   glyph=shortnames(ng.iris)[c(1,2,3,4,3,2,3,4,1)])
vizGlyph2


###################################################
### code chunk number 74: RnavGraph.Rnw:1070-1071 (eval = FALSE)
###################################################
## nav <- navGraph(ng.iris, list(ng.lg,ng.lgnot), list(vizGlyph1, vizGlyph2))


###################################################
### code chunk number 75: RnavGraph.Rnw:1111-1114
###################################################
myPlot <<- function(x,y,group,labels) {
  plot(x,y,col = group, pch = 19)
}


###################################################
### code chunk number 76: RnavGraph.Rnw:1118-1120
###################################################
viz1 <- ng_2d_myplot(ng.iris,ng.lg,fnName = "myPlot", device = "base")
viz1


###################################################
### code chunk number 77: RnavGraph.Rnw:1124-1125 (eval = FALSE)
###################################################
## nav <- navGraph(ng.iris,ng.lg, viz1)


###################################################
### code chunk number 78: RnavGraph.Rnw:1128-1129 (eval = FALSE)
###################################################
## dev.control(displaylist = "inhibit")


###################################################
### code chunk number 79: RnavGraph.Rnw:1142-1150 (eval = FALSE)
###################################################
## setClass(
##          Class="testVizClass",
##          representation =
##            representation(
##                           ## No additional slots
##                           ),
##          contains = "NG_Visualization"
##          )


###################################################
### code chunk number 80: RnavGraph.Rnw:1153-1167 (eval = FALSE)
###################################################
## myViz <- function(data,graph) {
##   if(is(data,"NG_data") == FALSE){
##     stop("data is no NG_data object.\n")
##   }
##   if(is(graph,"NG_graph") == FALSE){
##     stop("graph is no NG_graph object.\n")
##   }
##   
##   return(new(
##              "testVizClass",
##              graph = graph@name,
##              data = data@name
##              ))
## }


###################################################
### code chunk number 81: RnavGraph.Rnw:1170-1177 (eval = FALSE)
###################################################
## setMethod(
##           f = "initializeViz",
##           signature = "testVizClass",
##           definition = function(viz,ngEnv){
##             print(paste('You switched to the graph', viz@graph))
##             return(viz)
##           })


###################################################
### code chunk number 82: RnavGraph.Rnw:1180-1190 (eval = FALSE)
###################################################
## setMethod(
##           f = "updateViz",
##           signature = "testVizClass",
##           definition = function(viz,ngEnv){
##             print(paste('You current state is:', 
##                         ngEnv$bulletState$from, 'to',
##                         ngEnv$bulletState$to, 'and',
##                         floor(ngEnv$bulletState$percentage*100), 'percent in between'))
##             return(viz)
##           })


###################################################
### code chunk number 83: RnavGraph.Rnw:1193-1200 (eval = FALSE)
###################################################
## setMethod(
##           f = "closeViz",
##           signature = "testVizClass",
##           definition = function(viz,ngEnv){
##             print(paste('Bye Bye', viz@graph))
##             return(viz)
##           })


###################################################
### code chunk number 84: RnavGraph.Rnw:1204-1205 (eval = FALSE)
###################################################
## demo("ng_own_viz")


###################################################
### code chunk number 85: RnavGraph.Rnw:1334-1336 (eval = FALSE)
###################################################
## navGraph(..., settings=list(color=list(background="green"),
##                 interaction=list(bulletRadius=4, nodeRadius=3)))


###################################################
### code chunk number 86: RnavGraph.Rnw:1430-1433 (eval = FALSE)
###################################################
## nav <- navGraph(ng.iris)
## 
## ng_walk(nav,"s.L:s.W s.L:p.L p.L:p.W s.L:p.W")


###################################################
### code chunk number 87: RnavGraph.Rnw:1436-1437 (eval = FALSE)
###################################################
## ng_walk(nav,c("s.L:s.W","s.L:p.L","p.L:p.W","s.L:p.W"))


###################################################
### code chunk number 88: RnavGraph.Rnw:1445-1456
###################################################
library(RnavGraph)

## import the data
ng.iris <- ng_data(name = "iris", data = iris[,1:4],
		shortnames = c('s.L', 's.W', 'p.L', 'p.W'),
		group = iris$Species,
		labels = substr(iris$Species,1,2))

## start navGraph
nav <- navGraph(ng.iris)
nav


###################################################
### code chunk number 89: RnavGraph.Rnw:1460-1462
###################################################
nav <- ng_update(nav)
nav


###################################################
### code chunk number 90: RnavGraph.Rnw:1465-1466
###################################################
ng_get(nav)


###################################################
### code chunk number 91: RnavGraph.Rnw:1469-1471
###################################################
ng_get(nav,"data")
ng_get(ng_get(nav,"data"),"group")[1:5]


###################################################
### code chunk number 92: RnavGraph.Rnw:1474-1475 (eval = FALSE)
###################################################
## nav1 <- navGraph(nav)


###################################################
### code chunk number 93: RnavGraph.Rnw:1481-1483 (eval = FALSE)
###################################################
## ng_walk(nav,"s.L:s.W s.L:p.L p.L:p.W s.L:p.W")
## ng_walk(nav,c("s.L:s.W","s.L:p.L","p.L:p.W","s.L:p.W"))


###################################################
### code chunk number 94: RnavGraph.Rnw:1528-1552
###################################################
data(olive)
ng.olive <- ng_data(name = "Olive",
		data = data.frame(olive[,-c(1,2)]))

nav <- navGraph(list(ng.iris,ng.olive)) ## Start navGraph Session with two data sets

## get colors
colIris <- ng_get_color(nav, "iris")
sizeOlive <- ng_get_size(nav, "Olive")

## the following will tell you to specify your data set
col <- ng_get_color(nav)

## set all points for iris data to red
ng_set_color(nav,"iris") <- "red"

## set to random colors
ng_set_color(nav,"iris") <- sample(c("red","yellow","green"),replace = TRUE, 150)

## set all points to same size
ng_set_size(nav,"iris") <- 1

## set set points to random sizes
ng_set_size(nav,'iris') <- sample(1:7, replace=TRUE, 150)


###################################################
### code chunk number 95: RnavGraph.Rnw:1558-1565
###################################################
nav <- ng_update(nav)
## For iris
attrIris <- ng_get(ng_get(nav,"data")[[1]],"group")
attrIris[1:4]
## For olive
attrOlive <- ng_get(ng_get(nav,"data")[[2]],"group")
attrOlive[1:4]


###################################################
### code chunk number 96: RnavGraph.Rnw:1572-1587 (eval = FALSE)
###################################################
## V <- shortnames(ng.iris)
## G <- completegraph(V)
## LG <- linegraph(G)
## LGnot <- complement(LG)
## 
## ng.lg <- ng_graph(name = '3D Transition', graph = LG, layout = 'circle')
## ng.lgnot <- ng_graph(name = '4D Transition', graph = LGnot, layout = 'circle')
## 
## viz1 <- ng_2d(ng.iris,ng.lg, glyphs = V[c(1,2,3,4,1,3,2,4)])
## viz2 <- ng_2d(ng.iris,ng.lg)
## 
## viz <- list(viz1, viz2)	
## graphs <- list(ng.lg, ng.lgnot)
## 
## nav <- navGraph(graph = graphs, data = ng.iris, viz = viz)


###################################################
### code chunk number 97: RnavGraph.Rnw:1593-1595 (eval = FALSE)
###################################################
## nav1 <- navGraph(ng.iris)
## nav2 <- navGraph(ng.iris)


###################################################
### code chunk number 98: RnavGraph.Rnw:1599-1600 (eval = FALSE)
###################################################
## nav3 <- navGraph(ng.iris, settings = list(tk2d=list(linked = FALSE)))


###################################################
### code chunk number 99: RnavGraph.Rnw:1608-1615
###################################################
library(scagnostics)
data(olive)
ng.olive <- ng_data(name = "Olive",
                    data = olive[,-c(1,2)],
                    shortnames = c("p1","p2","s","ol","l1","l2","a","e"),
                    group = as.numeric(olive$Area)+1
                    )


###################################################
### code chunk number 100: RnavGraph.Rnw:1620-1622
###################################################
scagMat <- scagnostics(olive)
rownames(scagMat)


###################################################
### code chunk number 101: RnavGraph.Rnw:1627-1633 (eval = FALSE)
###################################################
## nav <- scagNav(data = ng.olive,
##                scags = c("Skinny", "Sparse", "NotConvex"),
##                topFrac = 0.2,
##                combineFn = max,
##                glyphs = shortnames(ng.olive)[1:8],
##                sep = ':')


###################################################
### code chunk number 102: RnavGraph.Rnw:1638-1643 (eval = FALSE)
###################################################
## nav <- scagNav(data = ng.olive,
##                scags = c("Skinny", "Sparse", "NotConvex"),
##                topFrac = 0.2,
##                glyphs = shortnames(ng.olive)[1:8],
##                sep = '+')


###################################################
### code chunk number 103: RnavGraph.Rnw:1653-1656
###################################################
edgeWts <- scagEdgeWeights(data = ng.olive, scags = c("Clumpy", "Skinny"))
edgeWts$fromToEdgeMatrix[1:3,]
edgeWts$nodeNames


###################################################
### code chunk number 104: RnavGraph.Rnw:1659-1663
###################################################
edgeWts <- scagEdgeWeights(data = ng.olive,
                           scags = c("Clumpy", "Skinny"),
                           combineFn = max)
edgeWts$fromToEdgeMatrix[1:3,]


###################################################
### code chunk number 105: RnavGraph.Rnw:1666-1670
###################################################
edgeWts <- scagEdgeWeights(data = ng.olive,
                           scags = c("Clumpy", "Skinny"),
                           combineFn = function(x){2*x[1]+3*x[2]})
edgeWts$fromToEdgeMatrix[1:3,]


###################################################
### code chunk number 106: RnavGraph.Rnw:1678-1683
###################################################
weights <- edgeWts$fromToEdgeMatrix[,"combined weights"]
ii <- weights>quantile(weights,0.8)
G <- newgraph(nodeNames = edgeWts$nodeNames,
		mat = edgeWts$fromToEdgeMatrix[ii,c(1,2)], weights = weights[ii])
plot(G)


###################################################
### code chunk number 107: RnavGraph.Rnw:1688-1689
###################################################
edgeData(G, attr="weight")$'p1|e'


###################################################
### code chunk number 108: RnavGraph.Rnw:1693-1695
###################################################
ng.lg <- ng_graph("3d olive",linegraph(G))
viz <- ng_2d(ng.olive,ng.lg)


###################################################
### code chunk number 109: RnavGraph.Rnw:1697-1698 (eval = FALSE)
###################################################
## nav <- navGraph(ng.olive,ng.lg,viz)


###################################################
### code chunk number 110: RnavGraph.Rnw:1705-1715
###################################################
par(mfrow = c(2,2))
G_1 <- scagGraph(edgeWts, topFrac = 0.2)
plot(G_1)
edgeData(G_1,attr = "weight")$'p1|e'
G_1 <- scagGraph(edgeWts, topFrac = 0)
plot(G_1)
G_1 <- scagGraph(edgeWts, topFrac = 1)
plot(G_1)
G_1 <- scagGraph(edgeWts, topFrac = 0.0001)
plot(G_1)


###################################################
### code chunk number 111: RnavGraph.Rnw:1719-1725
###################################################
edgeWts <- scagEdgeWeights(data = ng.olive,
		scags = c("Clumpy","NotClumpy","Monotonic"),		
		combineFn = NULL)

graphList <- scagGraph(edgeWts, topFrac = 0.2)
graphList


###################################################
### code chunk number 112: RnavGraph.Rnw:1729-1730 (eval = FALSE)
###################################################
## nav <- navGraph(ng.olive, graphList)


###################################################
### code chunk number 113: RnavGraph.Rnw:1737-1738 (eval = FALSE)
###################################################
## library(MASS)


###################################################
### code chunk number 114: RnavGraph.Rnw:1742-1748 (eval = FALSE)
###################################################
## ng.iris <- ng_data(name = "iris", data = iris[,1:4],
## 		shortnames = c('s.L', 's.W', 'p.L', 'p.W'),
## 		group = as.numeric(iris$Species),
## 		labels = substr(iris$Species,1,2))
## 
## nav <- navGraph(ng.iris)


###################################################
### code chunk number 115: RnavGraph.Rnw:1752-1774 (eval = FALSE)
###################################################
## library(PairViz)
## data(olive)
## 
## d.olive <- data.frame(olive[,-c(1,2)])
## ng.olive <- ng_data(name = "Olive",
## 		data = d.olive,
## 		shortnames = c("p1","p2","s","oleic","l1","l2","a","e"),
## 		group = as.numeric(olive[,"Area"]),
## 		labels = as.character(olive[,"Area"])
## )
## ng.olive
## 
## G <- completegraph(shortnames(ng.olive))
## LG <- linegraph(G)
## ng.lg <- ng_graph("3d olive",LG, layout = 'kamadaKawaiSpring')
## ng.lgnot <- ng_graph("4d olive",complement(LG), layout = 'kamadaKawaiSpring')
## 
## nav <- navGraph(ng.olive,
## 		list(ng.lg,ng.lgnot),
## 		list(ng_2d(ng.olive,ng.lg),ng_2d(ng.olive,ng.lgnot)))
## 
## ng_walk(nav, eulerian(as(LG,"graphNEL"))[1:7])  


###################################################
### code chunk number 116: RnavGraph.Rnw:1778-1790 (eval = FALSE)
###################################################
## library(MASS)
## ng.data <- ng_data(name = "US Judge Ratings", data = USJudgeRatings)
## 
## p <- ncol(USJudgeRatings)
## adjM <- matrix(0,ncol=p,nrow=p)
## adjM[c(1:8,11),c(9,10,12)] <- 1
## adjM[c(9,10,12),c(1:8,11)] <- 1
## 
## G <- newgraph(names(ng.data),adjM, isAdjacency = TRUE)
## ng.lg <- ng_graph("3d Us Judge",linegraph(G),layout = 'fruchtermanReingold')
## 
## nav <- navGraph(ng.data,ng.lg,ng_2d(ng.data,ng.lg))


###################################################
### code chunk number 117: RnavGraph.Rnw:1795-1821 (eval = FALSE)
###################################################
## names(stormtracks)
## storms <- stormtracks[,c(2:9,11)]
## ng.storms <- ng_data(name = "Storm tracks",
## 		data = stormtracks[,c(2:9,11)],
## 		group = as.numeric(stormtracks[,"type"]),
## 		labels = stormtracks[,"type"])
## 
## p <- ncol(ng_get(ng.storms,"data"))
## adjM <- matrix(0,ncol=p,nrow=p)
## adjM[c(1,2,4,5,6),c(7,8,9)] <- 1
## adjM[c(7,8,9),c(1,2,4,5,6)] <- 1
## adjM[c(7,8,9),c(7,8,9)] <- 1
## adjM[7,7] <- adjM[8,8] <- adjM[9,9] <- 0
## adjM[c(5,6),c(5,6)] <- 1
## adjM[5,5] <- adjM[6,6] <- 0
## 
## G <- newgraph(names(ng.storms),adjM, isAdjacency = TRUE)
## LG <- linegraph(G)
## ng.lg <- ng_graph("3d storm tracks", LG, layout = 'kamadaKawaiSpring')
## ng.lgnot <- ng_graph("4d storm tracks", complement(LG), layout = 'kamadaKawaiSpring')
## viz1 <- ng_2d(ng.storms,ng.lg)
## viz2 <- ng_2d(ng.storms,ng.lgnot)
## 
## nav <- navGraph(ng.storms,
## 		list(ng.lg,ng.lgnot),
## 		list(viz1,viz2))


###################################################
### code chunk number 118: RnavGraph.Rnw:1826-1834 (eval = FALSE)
###################################################
## ng.data <- ng_data(name = "USCereal",
## 		data = UScereal[,c(2:8,10)],
## 		shortnames = c("cal","prot","fat","sod","fib","carb","sug","pt"),
## 		group = UScereal[,1],
## 		labels = UScereal[,1])
## 
## nav <- navGraph(ng.data)
## nav <- scagNav(ng.data, scags = "Outlying", topFrac = 0.2)


###################################################
### code chunk number 119: RnavGraph.Rnw:1838-1845 (eval = FALSE)
###################################################
## ng.data <- ng_data(name = "Boston",
## 		data = Boston[,-9],
## 		shortnames = names(Boston[,-9]),
## 		group = Boston[,"rad"])
## 
## nav <- navGraph(ng.data)
## nav <- scagNav(ng.data, scags = "Clumpy", topFrac = 0.2)


###################################################
### code chunk number 120: RnavGraph.Rnw:1849-1864 (eval = FALSE)
###################################################
## ng.data <- ng_data(name = "Birth Weight Data",
## 		data = birthwt[,c(1:3,5:10)],
## 		group = birthwt[,4])
## 
## p <- ncol(ng_get(ng.data,"data"))
## adjM <- matrix(0,ncol=p,nrow=p)
## adjM[c(1:8),c(9)] <- 1
## adjM[c(9),c(1:8)] <- 1
## 
## G <- newgraph(names(ng.data),adjM,isAdjacency = TRUE)
## LG <- linegraph(G, sep = "++")
## ng.lg <- ng_graph("3d birth weight", LG, sep = "++")
## ng.lgnot <- ng_graph("34 birth weight", complement(LG), sep = "++")
## 
## nav <- navGraph(ng.data, list(ng.lg, ng.lgnot), list(ng_2d(ng.data,ng.lg),ng_2d(ng.data,ng.lgnot)))


###################################################
### code chunk number 121: RnavGraph.Rnw:1868-1878 (eval = FALSE)
###################################################
## require(alr3)
## data(banknote)
## names(banknote[,1:6])
## ng.data <- ng_data(name = "Swiss bank note Data",
## 		data = banknote[,1:6],
## 		shortnames = names(banknote[,1:6]),
## 		group = banknote[,7])
## 
## nav <- navGraph(ng.data)
## nav <- scagNav(ng.data, scags = "Clumpy", topFrac = 0.2)


###################################################
### code chunk number 122: RnavGraph.Rnw:1882-1891 (eval = FALSE)
###################################################
## require(gclus)
## data(body)
## names(body[,1:10])
## ng.data <- ng_data(name = "Body Dimensions",
## 		data = body[,1:24],
## 		group = body[,25])
## 
## nav <- navGraph(ng.data)
## nav <- scagNav(ng.data, scags = "Clumpy", topFrac = 0.1)


###################################################
### code chunk number 123: RnavGraph.Rnw:1895-1914 (eval = FALSE)
###################################################
## require(gclus)
## data(ozone)
## ng.data <- ng_data(name = "Ozone data",
## 		data = ozone)
## 
## nav <- navGraph(ng.data)
## 
## ## Breiman and Friedman's selection
## p <- ncol(ng_get(ng.data,"data"))
## adjM <- matrix(0,ncol=p,nrow=p)
## adjM[c(1),c(2,4,5,6)] <- 1
## adjM[c(2,4,5,6),c(1)] <- 1
## 
## G <- newgraph(names(ng.data),adjM, isAdjacency = TRUE)
## LG <- linegraph(G)
## ng.lg <- ng_graph("3d ozone", LG, layout="circle")
## ng.lgnot <- ng_graph("4d ozone", complement(LG), layout="circle")
## 
## nav <- navGraph(ng.data, list(ng.lg, ng.lgnot), list(ng_2d(ng.data,ng.lg), ng_2d(ng.data,ng.lgnot)))


###################################################
### code chunk number 124: RnavGraph.Rnw:1918-1933 (eval = FALSE)
###################################################
## ng.data <- ng_data(name = "SwissFertility",
## 		data = swiss,
## 		shortnames = c("Fer","Agri","Exam","Edu","Cath","IM"))
## 
## p <- ncol(swiss)
## adjM <- matrix(0,ncol=p,nrow=p)
## adjM[1:5,6] <- 1
## adjM[6,1:5] <- 1
## 
## G <- newgraph(shortnames(ng.data),adjM,isAdjacency=TRUE)
## LG <- linegraph(G)
## ng.lg <- ng_graph("3d fertility", LG, layout='circle')
## ng.lgnot <- ng_graph("4d fertility", complement(LG), layout="fruchtermanReingold")
## 
## nav <- navGraph(ng.data, list(ng.lg, ng.lgnot), list(ng_2d(ng.data,ng.lg), ng_2d(ng.data,ng.lgnot)))


###################################################
### code chunk number 125: RnavGraph.Rnw:1937-1956 (eval = FALSE)
###################################################
## require(alr3)
## data(challeng)
## 
## ng.data <- ng_data(name = "Challenger Data",
## 		data = challeng[,c(1:3,5:7)])
## 
## 
## p <- ncol(ng_get(ng.data,"data"))
## adjM <- matrix(0,ncol=p,nrow=p)
## adjM[c(1,2),c(3:6)] <- 1
## adjM[c(3:6),c(1,2)] <- 1
## 
## G <- newgraph(names(ng.data),adjM,isAdjacency = TRUE)
## LG <- linegraph(G)
## 
## ng.lg <- ng_graph("3d challenger", LG)
## ng.lgnot <- ng_graph("4d challenger", complement(LG))
## 
## nav <- navGraph(ng.data, list(ng.lg, ng.lgnot), list(ng_2d(ng.data,ng.lg), ng_2d(ng.data,ng.lgnot)))


###################################################
### code chunk number 126: RnavGraph.Rnw:1960-1970 (eval = FALSE)
###################################################
## library(PairViz)
## require(cluster)
## data(animals)
## names(animals)
## ng.data <- ng_data(name = "Animal Data",
## 		data = animals)
## 
## nav <- navGraph(ng.data)
## 
## ng_walk(nav, eulerian(as(ng_get(ng_get(nav,"graphs")[[1]],"graph")),"graphNEL"))


