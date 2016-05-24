## ----setup, echo=FALSE---------------------------------------------------
# set global chunk options: images will be 7x5 inches
knitr::opts_chunk$set(fig.width=7, fig.height=7, fig.path="figs/")
options(digits = 4)

## ----install, eval=FALSE-------------------------------------------------
#  library(devtools)
#  install_github("thibautjombart/treescape")

## ----install2, eval=FALSE------------------------------------------------
#  install.packages("treescape")

## ----load----------------------------------------------------------------
library("treescape")

## ------------------------------------------------------------------------
library("treescape")
library("ade4")
library("adegenet")
library("adegraphics")
library("ggplot2")

## ----treescape-----------------------------------------------------------
## generate list of trees
set.seed(1)
x <- rmtree(10, 20)
names(x) <- paste("tree", 1:10, sep = "")

## use treescape
res <- treescape(x, nf=3)
names(res)
res

## ----distances-----------------------------------------------------------
## table.image
table.image(res$D, nclass=30)

## table.value with some customization
table.value(res$D, nclass=5, method="color", 
            symbol="circle", col=redpal(5))


## ----treescapescatter----------------------------------------------------
scatter(res$pco)

## ----plotgroves----------------------------------------------------------
plotGroves(res$pco, lab.show=TRUE, lab.cex=1.5)

## ----woodmicePlots-------------------------------------------------------
data(woodmiceTrees)
wm.res <- treescape(woodmiceTrees,nf=3)

## this is the PCoA / MDS:
wm.res$pco

## PCs are stored in:
head(wm.res$pco$li)

## plot results
plotGroves(wm.res$pco, lab.show=TRUE, lab.optim=FALSE)

## visualising density of points
s.kde2d(wm.res$pco$li)

## alternative visualisation
s.density(wm.res$pco$li, col=redpal(100), bandwidth=3)

## same, other palette
s.density(wm.res$pco$li, col=rev(transp(spectral(100),.5)), bandwidth=3)

## alternative using ggplot2
woodmiceplot <- ggplot(wm.res$pco$li, aes(x=A1, y=A2)) # create plot
woodmiceplot + geom_density2d(colour="gray80") + # contour lines
geom_point(size=6, shape=1, colour="gray50") + # grey edges
geom_point(size=6, alpha=0.2, colour="navy") + # transparent blue points
xlab("") + ylab("") + theme_bw(base_family="") # remove axis labels and grey background

## ----findgroves, cache=TRUE----------------------------------------------
wm.groves <- findGroves(wm.res, nclust=6)
names(wm.groves)

## ----plotgroves2---------------------------------------------------------
## basic plot
plotGroves(wm.groves)

## alternative with inertia ellipses
plotGroves(wm.groves, type="ellipse")

## plot axes 2-3
plotGroves(wm.groves, xax=2, yax=3)

## ----plotgroves3, fig.keep="last"----------------------------------------
## customize graphics
plotGroves(wm.groves, bg="black", col.pal=lightseasun, lab.show=TRUE, lab.col="white", lab.cex=1.5)


## ----woodmiceMedian------------------------------------------------------
## get first median tree
tre <- medTree(woodmiceTrees)$trees[[1]]

## plot tree
plot(tre,type="cladogram",edge.width=3, cex=0.8)

## ----woodmiceCluster1, out.width="600px"---------------------------------
## find median trees for the 6 clusters identified earlier:
res <- medTree(woodmiceTrees, wm.groves$groups)

## there is one output per cluster
names(res)

## get the first median of each
med.trees <- lapply(res, function(e) ladderize(e$trees[[1]]))

## plot trees
par(mfrow=c(2,3))
for(i in 1:length(med.trees)) plot(med.trees[[i]], main=paste("cluster",i),cex=1.5)


## ----woodmice-tip-emphasis-----------------------------------------------
wm3.res <- treescape(woodmiceTrees,nf=2,emphasise.tips=c("No1007S","No1208S","No0909S"),emphasise.weight=3)

## plot results
plotGroves(wm3.res$pco, lab.show=TRUE, lab.optim=FALSE)

## ----findgroves-with-emphasis--------------------------------------------
wm3.groves <- findGroves(woodmiceTrees,nf=3,nclust=6,emphasise.tips=c("No1007S","No1208S","No0909S"),emphasise.weight=3)
plotGroves(wm3.groves, type="ellipse")

## ----treevec-------------------------------------------------------------
## generate a random tree:
tree <- rtree(6)
## topological vector of mrca distances from root:
treeVec(tree)
## vector of mrca distances from root when lambda=0.5:
treeVec(tree,0.5)
## vector of mrca distances as a function of lambda:
vecAsFunction <- treeVec(tree,return.lambda.function=TRUE)
## evaluate the vector at lambda=0.5:
vecAsFunction(0.5)

## ----treedist------------------------------------------------------------
## generate random trees
tree_a <- rtree(6)
tree_b <- rtree(6)

## topological (lambda=0) distance:
treeDist(tree_a,tree_b) 

## branch-length focused (lambda=1) distance:
treeDist(tree_a,tree_b,1)

