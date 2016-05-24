### R code from vignette source 'catnet.Rnw'

###################################################
### code chunk number 1: cnet
###################################################
library(catnet)
cnet <- cnNew(
	nodes = c("a", "b", "c"),
	cats = list(c("1","2"), c("1","2"), c("1","2")), 
	parents = list(NULL, c(1), c(1,2)), 
	probs = list(	c(0.2,0.8), 
		list(c(0.6,0.4),c(0.4,0.6)), 
		list(list(c(0.3,0.7),c(0.7,0.3)), 
		list(c(0.9,0.1),c(0.1,0.9))))  )


###################################################
### code chunk number 2: cnet1
###################################################
set.seed(123)
cnet1 <- cnRandomCatnet(numnodes=4, maxParents=2, numCategories=2)
cnet1


###################################################
### code chunk number 3: cnet2
###################################################
myNodes<-c("a","s","p","q","r","t","u")
myEdges<-list(a=list(edges=NULL), s=list(edges=c("p","q")), p=list(edges=c("q")), q=list(edges=c("r")), r=list(edges=c("u")), t=list(edges=c("q")),u=list(edges=NULL))
cnet2 <- cnCatnetFromEdges(nodes=myNodes, edges=myEdges, numCategories=2)


###################################################
### code chunk number 4: cnet3
###################################################
#cnet3 <- cnCatnetFromSif(filename)
#cnet3


###################################################
### code chunk number 5: simplefuns
###################################################
cnNumNodes(cnet1)
cnNodes(cnet1)
cnEdges(cnet1)
cnParents(cnet1)


###################################################
### code chunk number 6: simplefuns2
###################################################
cnMatParents(cnet1)
cnMatEdges(cnet1)


###################################################
### code chunk number 7: simplefuns3
###################################################
cnPlotProb(cnet1)


###################################################
### code chunk number 8: simplefuns4
###################################################
cnComplexity(cnet1)


###################################################
### code chunk number 9: ex2
###################################################
cnOrder(cnet1)
cnOrder(cnet1@parents)


###################################################
### code chunk number 10: ex4
###################################################
set.seed(456)
cnet2 <- cnRandomCatnet(numnodes=10, maxParents=3, numCategories=2)
cnEdges(cnet2)
pcnet2 <- dag2cpdag(cnet2)
cnEdges(pcnet2)


###################################################
### code chunk number 11: ex5
###################################################
set.seed(456)
cnet3 <- cnRandomCatnet(cnNumNodes(cnet2), maxParents=2, numCategories=2)
cnet3@nodes <- cnet2@nodes
cnCompare(object1=cnet2, object2=cnet3)


###################################################
### code chunk number 12: ex6
###################################################
samples1 <- cnSamples(object=cnet1, numsamples=100, output="matrix")
dim(samples1)
samples1 <- cnSamples(object=cnet1, numsamples=100, output="frame")
dim(samples1)


###################################################
### code chunk number 13: ex7
###################################################
samples2 <- cnSamples(object=cnet1, numsamples=10, perturbations=c(0,0,1,2))


###################################################
### code chunk number 14: ex8
###################################################
  ## generate a sample of size 12 and set the last 3 nodes as not-available
  numnodes <- cnNumNodes(cnet2)
  samples3 <- cnSamples(object=cnet2, numsamples=12, output="matrix")
  ## predict the last three nodes in 'cnet2' from the rest
  ## by setting their values in 'samples3' as NA
  samples3[numnodes-2, ] <- rep(NA, 12)
  samples3[numnodes-1, ] <- rep(NA, 12)
  samples3[numnodes, ] <- rep(NA, 12)
  ## predict the values of the last 3 nodes
  newsamples <- cnPredict(object=cnet2, data=samples3)


###################################################
### code chunk number 15: ex12
###################################################
set.seed(789)
cnet2 <- cnRandomCatnet(numnodes=10, maxParents=2, numCategories=2)
nodeOrder <- order(runif(cnNumNodes(cnet2)))
cnet2
## generate a 100-size sample from cnet2
samples <- cnSamples(object=cnet2, numsamples=100, output="frame")
netlist <- cnSearchOrder(data=samples, perturbations=NULL, maxParentSet=2, maxComplexity=20, 
	nodeOrder, parentsPool=NULL, fixedParents=NULL)
## find the recostructed network with the true complexity
bnet <- cnFind(netlist, 20)
bnet


###################################################
### code chunk number 16: ex14
###################################################
set.seed(123)
nnodes <- 12
cnet <- cnRandomCatnet(numnodes=nnodes, maxParents=5, numCategories=2)
norder <- cnOrder(cnet)
parPool <- vector("list", nnodes)
for(i in 1:nnodes) parPool[[i]] <- 1:(nnodes-1)
fixparPool <- vector("list", nnodes)
for(i in 3:nnodes) fixparPool[[i]] <- c(1,2)
samples <- cnSamples(cnet, numsamples=200)
eval <- cnSearchOrder(data=samples, perturbations=NULL, maxParentSet=2, maxComplexity=200, 
	nodeOrder=norder, parentsPool=parPool, fixedParents=fixparPool)
eval
## plot likelihood vs complexity for the resulting list of networks
## plot(eval@complexity, eval@loglik, xlab="Complexity", ylab = "Log-likelihood", main="Model selection curve for the list of networks")


###################################################
### code chunk number 17: ex20
###################################################
set.seed(345)
## generate a 100-size sample from cnet6
cnet6 <- cnRandomCatnet(numnodes=12, maxParents=5, numCategories=2)
samples <- cnSamples(object=cnet6, numsamples=100, output="matrix")
eval <- cnSearchOrder(data=samples, perturbations=NULL, maxParentSet=2, parentSizes=NULL, maxComplexity=0, 
	nodeOrder=order(runif(1:dim(samples)[1])), parentsPool=NULL, fixedParents=NULL, 
	echo=FALSE)
## now select a network based on AIC and plot it
anet <- cnFindAIC(object=eval)
anet
## or BIC
bnet <- cnFindBIC(object=eval, numsamples=dim(samples)[2])
bnet
## plot likelihood vs complexity for the resulting list of networks
plot(eval@complexity, eval@loglik, 
xlab="Complexity", ylab = "Log-likelihood", 
main="Model selection: AIC and BIC complexities in red and blue.")
abline(v=anet@complexity,lty=2,col="red")
abline(v=bnet@complexity,lty=3,col="blue")


###################################################
### code chunk number 18: ex18
###################################################
set.seed(345)
samples <- cnSamples(object=cnet6, numsamples=100, output="matrix")
netlist <- cnSearchSA(data=samples, perturbations=NULL, 
	maxParentSet=2, parentSizes=NULL, maxComplexity=20, 
	parentsPool=NULL, fixedParents=NULL, 
	tempStart=1, tempCoolFact = 0.9, tempCheckOrders = 4, maxIter = 40, 
	orderShuffles = 1, stopDiff = 0.0001,
	priorSearch=NULL)
bnet <- cnFind(netlist@nets, cnComplexity(cnet6))
bnet


###################################################
### code chunk number 19: ex20
###################################################
set.seed(678)
numnodes <- 16
numcats <- 3
maxpars <- 2
cnet8 <- cnRandomCatnet(numnodes, maxpars, numcats)
ps <- cnSamples(cnet8, 500)
## next, a variable number of categories scanario is demonstrated
## find a node with descendants and reduce its number of categories
mpars <- cnMatParents(cnet8)
for(j in 1:numnodes)
	if(sum(mpars[,j]) > 0)
		break
if(j < numnodes)
	cnet8@categories[[j]] <- cnet8@categories[[j]][1:(numcats-1)]
## now resets cnet8's probability table
cnet8 <- cnSetProb(cnet8, ps)
## generate a new sample from the updated network
ps <- cnSamples(cnet8, 500)
res8 <- cnSearchOrder(data=ps, perturbations = NULL,
                      maxParentSet = maxpars, parentSizes = NULL, 
                      maxComplexity = 0,
                      nodeOrder = cnOrder(cnet8),  
                      parentsPool = NULL, fixedParents = NULL,
                      edgeProb = NULL, 
                      echo = FALSE)
anet8 <- cnFind(res8, cnComplexity(cnet8))
cnCompare(cnet8, anet8)
## perform a stochastic search with a prior that favors the true edges (given by mpars)
edgeHisto <- 0.5 + mpars / 4
res9 <- cnSearchSA(data=ps, perturbations=NULL,
                    maxParentSet=1, parentSizes = NULL, 
                    maxComplexity = 0,
                    parentsPool = NULL, fixedParents = NULL, edgeProb = edgeHisto, 
                    selectMode = "BIC", 
                    tempStart = 1, tempCoolFact = 0.9, tempCheckOrders = 20, 
                    maxIter = 100, orderShuffles = -1, stopDiff = 1,
                    numThreads = 2, 
                    priorSearch = NULL,
                    echo=FALSE)
anet9 <- cnFind(res9, cnComplexity(cnet8))
cnCompare(cnet8, anet9)


###################################################
### code chunk number 20: ex21
###################################################
cnSetSeed(3456)
ncats <- 2
cn <- cnRandomCatnet(20, 3, ncats, p.delta1=0.1, p.delta2=0.2) 
norder <- cnOrder(cn)
numnodes <- cnNumNodes(cn)
mpars <- cnMatParents(cn)
numsamples <- 100
## simulate perturbations
pert <- as.data.frame(matrix(rbinom(numnodes*numsamples, 1, p=0.25), ncol=numnodes))
for(j in 1:numsamples) 
  for(i in 1:numnodes) {
    if(pert[j,i])
      pert[j,i] <- 1+floor(runif(1)*ncats)
    for(ip in cn@parents[[i]]) {
      if(pert[j,ip]) {
        pert[j,i] <- 0
      }
    }
  }
ps <- cnSamples(cn, numsamples, pert, as.index=TRUE)
## pairwise conditional entropy difference between perturbed and non-perturbed data
klmat <- cnEdgeDistanceKL(ps, pert)
fscore1 <- NULL
fscore2 <- NULL
for(ntrials in 1:5) {
  numiter <- 60
  sares1 <- cnSearchSA(data=ps, perturbations=pert, maxParentSet=2,
                       parentSizes=NULL, maxComplexity=0,
                       parentsPool=NULL, fixedParents=NULL,
                       edgeProb=NULL, dirProb=NULL, selectMode = "AIC",
                       tempStart=0.1, tempCoolFact=0.9, tempCheckOrders=numiter,
                       maxIter=numiter, 
                       orderShuffles=0, stopDiff=0,
                       numThreads=2, priorSearch=NULL, echo=FALSE)
  cmp <- cnCompare(cn, cnFind(sares1, cnComplexity(cn)))
  fscore1 <- c(fscore1, cmp@fscore)
  sares2 <- cnSearchSA(data=ps, perturbations=pert, maxParentSet=2,
                       parentSizes=NULL, maxComplexity=0,
                       parentsPool=NULL, fixedParents=NULL,
                       edgeProb=NULL, dirProb=t(klmat), selectMode = "AIC",
                       tempStart=0.1, tempCoolFact=0.9, tempCheckOrders=numiter,
                       maxIter=numiter, 
                       orderShuffles=0, stopDiff=0,
                       numThreads=2, priorSearch=NULL, echo=FALSE)
  cmp <- cnCompare(cn, cnFind(sares2, cnComplexity(cn)))
  fscore2 <- c(fscore2, cmp@fscore)
}
pl <- list("No Prior"=fscore1, "KL-dist Prior"=fscore2)
boxplot(pl, ylab="F-score")


