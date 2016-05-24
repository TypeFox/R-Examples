library(pcalg)

set.seed(123)
nreps <- 100
res <- logical(nreps)
all.eff.true <- res
Rnd <- function(e) round(e, 14)## get 14 digits accuracy, as we use true (DAG, cov)
for (i in 1:nreps) {
  p <- 2 + rpois(1, lambda = 8) # ==>  p >= 2, E[p] = 10
  ## generate and draw random DAG :
  myDAG <- randomDAG(p, prob = 0.2)
  myCPDAG <- dag2cpdag(myDAG)
  mcov <- trueCov(myDAG)

  ## x != y  in {1,2,...p} ;
  xy <- sample.int(p, 2); x <- xy[1]; y <- xy[2]

  ## plot(myCPDAG)
  eff.true <- Rnd(causalEffect(myDAG, y, x))
  all.eff.true[i] <- eff.true
  ## cat("x=",x," y=",y," eff=",eff.true,"\n")

  eff.est <- Rnd(ida(x,y, mcov, myCPDAG, method="local"))
  res[i] <- (eff.true %in% eff.est)
}
cat('Time elapsed: ', (.pt <- proc.time()),"\n")

stem(all.eff.true)
if (!all(res)) stop("Test ida: True effects were not recovered!")

## *one* test for  method="global" :
eff.g.est <- Rnd(ida(x,y, mcov, myCPDAG, method="global", verbose=TRUE))
stopifnot(eff.est == eff.g.est)

cat('Time elapsed additionally: ', proc.time() - .pt,"\n")

## another special case (from Raphael Gervais)
set.seed(123)
p <- 7
myDAG <- randomDAG(p, prob = 0.2) ## true DAG
amatT <- as(myDAG, "matrix") # weighted adjacency matrix of true DAG
effT <- Rnd(amatT[2,3]*amatT[3,5]) # Causal effect of 2 on 5 from true DAG weighted matrix
myCPDAG <- dag2cpdag(myDAG) ## true CPDAG
covTrue <- trueCov(myDAG) ## true covariance matrix
effG <- Rnd(ida(2,5, covTrue,myCPDAG,method = "global"))

if (!(effT %in% effG)) stop("Test ida special case: True effects were not recovered!")
