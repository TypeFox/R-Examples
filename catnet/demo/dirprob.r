library(catnet)

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

## find the node conditional entropy difference between the non-perturbed and perturbed samples
klmat <- cnEdgeDistanceKL(ps, pert)

## stochastic search without and with directional priors 
fscore1 <- NULL
fscore2 <- NULL
for(ntrials in 1:10) {
  numiter <- 100
  sares1 <- cnSearchSA(data=ps, perturbations=pert, maxParentSet=3,
                       parentSizes=NULL, maxComplexity=0,
                       parentsPool=NULL, fixedParents=NULL,
                       edgeProb=NULL, dirProb=NULL, selectMode = "AIC",
                       tempStart=0.1, tempCoolFact=0.9, tempCheckOrders=numiter,
                       maxIter=numiter, 
                       orderShuffles=0, stopDiff=0,
                       numThreads=2, priorSearch=NULL, echo=FALSE)
  cmp <- cnCompare(cn, cnFind(sares1, cnComplexity(cn)))
  fscore1 <- c(fscore1, cmp@fscore)

  sares2 <- cnSearchSA(data=ps, perturbations=pert, maxParentSet=3,
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

## compare the two experiments
pl <- list("No Prior"=fscore1, "KL-dist Prior"=fscore2)
boxplot(pl, ylab="F-score")
