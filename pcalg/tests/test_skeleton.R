library(pcalg)
source(system.file(package="Matrix", "test-tools-1.R", mustWork=TRUE))
##--> showProc.time(), assertError(), relErrV(), ...

set.seed(234)
p <- 10
nreps <- 30
resG <- resP <- rep(FALSE,nreps)
for (i in 1:nreps) {
  ## generate and draw random DAG :
  myDAG <- randomDAG(p, prob = 0.3)
  mycor <- cov2cor(trueCov(myDAG))
  amat <- wgtMatrix(myDAG)
  amat[amat!=0] <- 1
  amat <- amat + t(amat)
  amat[amat!=0] <- 1

  ## Gaussian
  suffStat <- list(C = cov2cor(trueCov(myDAG)), n = 10^9)
  indepTest <- gaussCItest

  resU <- skeleton(suffStat, indepTest, 0.99, p = p,
      method = ifelse(i < nreps/2, "stable", "stable.fast"))

  resG[i] <- all(as(resU@graph,"matrix") == amat)
  resP[i] <- all(resU@pMax[as(resU@graph,"matrix") == TRUE] < 0.99)
}
showProc.time()

if (!all(resG)) stop("Test skeleton wrong: Some skeleton was not found correctly!")
if (!all(resP)) stop("Test skeleton wrong: There was an inconsistency with an entry in pMax!")

(doExtras <- pcalg:::doExtras())
if(!doExtras && !interactive()) q("no")

## ../inst/xtraR/graph2ftmatrix.R :
source(system.file(package="pcalg", "xtraR", "graph2ftmatrix.R", mustWork=TRUE))
##-> graph2ftM(), perm.ftM()

##' checks the simulation results below, assuming that res[[1]] is *un*permuted
isEqPerms <- function(res) {
    stopifnot((M <- length(res)) > 1, vapply(res, is.list, NA),
              is.list(sk1 <- res[[1]]$skels))
    sapply(names(sk1), function(jj) {
        g1 <- graph2ftM(sk1[[jj]]@graph)
        vapply(2:M, function(m)
               identical(perm.ftM(g1, res[[m]]$perm),
                         graph2ftM(res[[m]]$skels[[jj]]@graph)), NA)
    })
}

showSkel <- function(skels) {
    stopifnot(is.list(skels))
    cat("#{edgetests} :\n")         #
    print(lapply(skels, slot, "n.edgetests"))
    ## pp <- sapply(skels, slot, "pMax")
    gs <- lapply(skels, slot, "graph")
    cat("edges  i --> j :\n")
    print(try(lapply(gs, graph2ftM))) # [revealed buglet in Matrix <= 1.1-0]
}


alphas <- c(.9,.5,.1,.05,.02, .01,.005, .002, .001)
names(alphas) <- formatC(alphas)

## for setup experiments
nrep <- 17
p <- 9
n.perm <- 12

## for real
nrep <- 7
n.perm <- 5
p <- 100

showProc.time()
pr0 <- 1/ 2^ceiling(log2(p))
for(skel.meth in c("stable", "original", "stable.fast")) {
    lin <-		"=============================="
    cat(sprintf("\n\n%s\nSkeleton method = \"%s\"\n%s\n",
		lin, skel.meth, lin))
    for(i in 1:nrep) {
	set.seed(i)		      # so we can easily get into one example
	myDAG <- randomDAG(p, prob = min(3/4, pr0 + (i-1)*pr0/4))
	cat(sprintf("\n\ni = %2d; #{edges} = %3d\n-----------\n", #
		    i, numEdges(myDAG)))
	C <- cov2cor(trueCov(myDAG))
	res <- vector("list", n.perm)
	for(k in 1:n.perm) {
	    ## randomly permute the nodes (apart from first case)
	    ii <- if(k == 1) 1:p else sample.int(p)
	    ##	*not* a large n [so alpha makes a difference
	    sStat <- list(C = C[ii,ii], n = 100)
	    sk.s <- lapply(alphas, function(AA)
			   skeleton(sStat, indepTest=gaussCItest, p=p, alpha = AA,
				    method = skel.meth))
	    showSkel(sk.s)
	    res[[k]] <- list(perm = ii, skels = sk.s)
	}
	showProc.time()
	fname <- paste0("test_skeleton_", substr(skel.meth, 1,4),"_", i,".rda")
	cat("Saving finished skeletons() to '",fname,"':\n", sep='')
	save(myDAG, res, file=fname)
        ## Now check that  res[[.]] are all the "same" for PC-stable,
        ## but not for "PC-orig" :
	cat("Checking skeletons() ")
	switch(skel.meth,
	       "stable" = stopifnot(isEqPerms(res)),
         "stable.fast" = stopifnot(isEqPerms(res)),
	       "original" = {
		   eq <- isEqPerms(res)
		   if(!all(eq)) {
		       cat("Skeletons of permuted variables are not all the same!\n")
		       print(eq)
		   }
	       }, stop("invalid 'skel.meth': ", skel.meth))
	cat("[Ok]\n"); showProc.time()
    }
}


if(FALSE) { ## experiments

## (load("test_skeleton_6.rda")); res6 <- res; D6 <- myDAG; str(res6, max=1)
## (load("test_skeleton_17.rda")); res17 <- res; D17 <- myDAG; str(res17, max=1)

graph2ftM(D6)
et6 <- lapply(seq_along(res6), function(j)
              lapply(lapply(res6[[j]]$skels, slot, "graph"), graph2ftM))
## Number of edges:
(ne.6 <- sapply(et6, vapply, nrow, 1))
## are the same for each of n.perm  permutations:
stopifnot(all(ne.6[,1] == ne.6))

(P3 <- as(res6[[3]]$perm,"pMatrix"))
C6 <- cov2cor(trueCov(D6))

graph2ftM(res6[[1]]$skels[["0.05"]]@graph)

graph2ftM(res6[[2]]$skels[["0.05"]]@graph)
        res6[[2]]$perm
invPerm(res6[[2]]$perm)

graph2ftM(res6[[3]]$skels[["0.05"]]@graph)
        res6[[3]]$perm
invPerm(res6[[3]]$perm)

ftM <- graph2ftM(res6[[1]]$skels[["0.05"]]@graph)
perm <- res6[[2]]$perm

p <- length(ip <- invPerm(perm))
m <- ftM; m[] <- ip[ftM]
M <- graph2ftM(T2graph(new("nsTMatrix", i=m[,1]-1L, j=m[,2]-1L, Dim=c(p,p))))
M
}
