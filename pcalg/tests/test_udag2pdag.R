library(pcalg)

.libPaths()
## acyclic graphs

nreps <- 30
p <- 8
n <- 1000

for(u2pd in c("rand", "retry", "relaxed")) {
    cat("\n u2pd =", u2pd, "\n ------------\n")
    cyc.res <- logical(nreps)
    for (i in 1:nreps) {
        set.seed(i)
        myDAG <- randomDAG(p, prob = 0.2)
        d.mat <- rmvDAG(n, myDAG, errDist = "normal")
        res <- suppressWarnings(pcAlgo(d.mat, alpha = 0.05, directed=TRUE, u2pd = u2pd))
        ##                      ------ directed;  u2pd = "rand" --> udag2pdag()
        res.A <- wgtMatrix(res@graph)
        res.A[res.A!=0] <- 1
        undir.A <- res.A + t(res.A)
        undir.A[undir.A==1] <- 0
        undir.A[undir.A==2] <- 1
        res.dir <- res.A - undir.A
        cyc.res[i] <- ggm::isAcyclic(res.dir)
    }
    if (!all(cyc.res)) stop("Test of pcAlgo(*, directed): Cyclic part in PDAG!")
} ## for(u2pd ...)
cat('Time elapsed: ', (.pt <- proc.time()),"\n")

## find collider correctly
set.seed(123)
myDAG <- randomDAG(3, prob = 0.5)
library(Matrix)
as(myDAG,"sparseMatrix")
d.mat <- rmvDAG(n, myDAG, errDist = "normal")
res <- pcAlgo(d.mat, alpha = 0.05, corMethod = "standard",directed=TRUE)
gEst <- wgtMatrix(res@graph)
gTrue <- rbind(0,0, c(1,1,0))
if(!all(gEst==gTrue)) stop("Test of udag2pdag: Problem finding a collider!")

cat('Time elapsed: ', proc.time() - .pt,'\n') # "stats"

