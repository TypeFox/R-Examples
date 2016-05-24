library(pcalg)

set.seed(42)
acyc <- edgeDistr <- wgtDistr <- 0

p <- 10
s <- 0.2
lB <- 0
uB <- 2
wgts <- NULL
nEdg <- 0
nreps <- 10
acyc <- rep(FALSE,nreps)

for (i in 1:nreps) {
  myDAG <- randomDAG(p,s,lB,uB)
  amat <- as(myDAG,"matrix")
  acyc[i] <- ggm::isAcyclic(amat)
  wgts <- c(wgts,amat[amat!=0])
  nEdg <- nEdg + graph::numEdges(myDAG)
}

## Check number of edges
nPossibleEdg <- nreps*p*(p-1)/2
edge.test <- binom.test(nEdg,nPossibleEdg,p=s,alternative="two.sided")
if (edge.test$p.value < 0.05) {
  stop("Test of randomDAG: Number of edges is not as theory predicts!")
}

## Check distribution of wgts
if (min(wgts) < lB || max(wgts) > uB) {
  stop("Test of randomDAG: Weights are not within specified ranges!")
}

## Check whether graphs are acyclic
if (!all(acyc))
  stop("Test of randomDAG: Graph is not acyclic!")

## Check special case
if (graph::numEdges(randomDAG(n=2,prob=1)) != 1) {
    stop("Test of randomDAG: Error with special case n=2, prob=1!")
}

cat('Time elapsed: ', proc.time(),"\n")
