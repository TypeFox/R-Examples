### R code from vignette source 'geneNetBPvignette.Rnw'

###################################################
### code chunk number 1: geneNetBPvignette.Rnw:54-57
###################################################
require(geneNetBP)
prettyVersion <- packageDescription("geneNetBP")$Version
prettyDate <- format(Sys.Date())


###################################################
### code chunk number 2: dataset
###################################################
data(mouse,package="geneNetBP")
head(mousegeno,n=3)
head(mousepheno,n=3)


###################################################
### code chunk number 3: dataset (eval = FALSE)
###################################################
## data(yeast,package="geneNetBP")
## head(yeastgeno,n=3)
## head(yeastpheno,n=3)


###################################################
### code chunk number 4: B1 (eval = FALSE)
###################################################
## network<-fit.gnbp(mousegeno,mousepheno,alpha = 0.1)
## ##convert the RHugin domain to a graph object
## BNgraph<-as.graph.RHuginDomain(network$gp)
## ##set node font size
## attrs<-list()
## attrs$node$fontsize<-30
## ## plot method for graph objects
## plot(BNgraph,attrs=attrs) 


###################################################
### code chunk number 5: L (eval = FALSE)
###################################################
##   ## Load the toy dataset
##   data(toy)
## ## Create a list of edges ("from (parent)", "to (child)")
## edgelist=list()
## edgelist[[1]]<-cbind("Q1","X1")
## edgelist[[2]]<-cbind("Q2","X1")
## edgelist[[3]]<-cbind("Q2","X2")
## edgelist[[4]]<-cbind("Q2","X4")
## edgelist[[5]]<-cbind("X1","X2")
## edgelist[[6]]<-cbind("Q3","X2")
## edgelist[[7]]<-cbind("Q3","X3")
## edgelist[[8]]<-cbind("X2","X5")
## edgelist[[9]]<-cbind("X2","X6")
## edgelist[[10]]<-cbind("X4","X6")


