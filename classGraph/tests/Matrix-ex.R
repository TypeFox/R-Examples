###################################################
### chunk number 1: preliminaries
###################################################

library(classGraph)

stopifnot(require("Matrix"),
          require("graph"))

allCl <- getClasses("package:Matrix")
## Those called "...Matrix" :
M.Cl <- grep("Matrix$",allCl, value = TRUE)

str(M3cl <- grep("^...Matrix$",M.Cl, value = TRUE))

M3cl <- M3cl[M3cl != "corMatrix"] # corMatrix not desired in following





###################################################
### chunk number 2: tree-Matrix
###################################################
trMatrix <- classTree("Matrix") ; trMatrix
## 2005-08: 58 nodes with 109 edges
## 2007-01: 74 nodes with 151 edges
par(lwd = 0.3)# no effect --- learn to use Rgraphviz!!
plotRag(mRagraph(trMatrix), subArgs=.optRagargs(adj = 0.5))


###################################################
### chunk number 3: sTree-M1
###################################################
allN <- nodes(trMatrix)
"%w/o%" <- function(x,y) x[!x %in% y] #--  x without y

hier1 <- paste(c("diagonal", "triangular", "symmetric", "general", "comp"),
               "Matrix", sep = '') # composites can be factorized:  ^^^^

## dropping 1st ``dim'' hierarchy -- much less edges:
(trM1 <- subGraph(allN %w/o% hier1, trMatrix))
plotRag(mRagraph(trM1), subArgs=.optRagargs(adj = 0.5))


###################################################
### chunk number 4: sTree-top12
###################################################
defMatrix <- getClassDef("Matrix")
## Distances of subclasses:
table(subDist <- sapply(defMatrix@subclasses, slot, "distance"))
sub12 <- defMatrix@subclasses[subDist <= 2] # but not unique!
trM_top12 <- subGraph(c("Matrix", unique(names(sub12))), trMatrix)
plotRag(mRagraph(trM_top12)) ## first and second level virtual classes


###################################################
### chunk number 5: sTree-sparse
###################################################
trSpMatrix <- classTree("sparseMatrix")
trSpMatrix
plotRag(mRagraph(trSpMatrix, "dot"),
        subArgs= .optRagargs(side = 3, adj = 1, line = 0),
        main = "'sparseMatrix' classes -- sub graph")


###################################################
### chunk number 6: sTree-sp-neato
###################################################
## now this *does* look kind of neat:
plotRag(mRagraph(trSpMatrix, "neato"),
        subArgs= .optRagargs(side = 3, adj = 1, line = 0),
        main = "'sparseMatrix' classes -- in \"neato\" layout")
## all three others do not


