library(pcalg)

showProc.time <- local({
    pct <- proc.time()
    function() { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- proc.time()
	cat('Time elapsed: ', (pct - ot)[1:3],'\n')
    }
})


##################################################
## Standard PC
##################################################
##library(RBGL)
nreps <- 10

p <- 10
set.seed(234)
for (ii in 1:nreps) {
  ## generate and draw random DAG :
  myDAG <- randomDAG(p, prob = 0.3)
  myCPDAG <- dag2cpdag(myDAG)
  suffStat <- list(C = cov2cor(trueCov(myDAG)), n = 10^9)
  res <- pc(suffStat, indepTest=gaussCItest, 0.99, p=p)
  if( shd(res, myCPDAG) != 0)
    stop("Test pc wrong: CPDAG ",ii," was not found correctly")
}

showProc.time()

##################################################
## Conservative PC
##################################################
##PC algorithm sample (compared with Tetrad)
##__________________________________________________________________________________
## Example 1
p <- 10
n <- 5000
set.seed(p*37+15673)
g <- randomDAG(p,2/(p-1))

## generate n samples of DAG using standard normal error distribution
##random data

(load(system.file("external", "test_conservative_pc_data1.rda", package = "pcalg")))
## set.seed(67*37)
## new.mat <- rmvnorm(n,mean=rep(0,p),sigma=trueCov(g))
## save(new.mat, file = "/u/kalisch/research/packages/pcalg/inst/external/test_conservative_pc_data1.rda")
## load(file = "/u/kalisch/research/packages/pcalg/inst/external/test_conservative_pc_data1.rda")

suffStat.data <- list(C=cor(new.mat1),n=n)

##pcAlgo conservative sample
dag1 <- pc(suffStat.data, gaussCItest, alpha=0.005, p=p,
	   u2pd="relaxed", conservative=TRUE)
##adjacency matrix
dag1.amat <- as(dag1@graph,"matrix")

##always save the transpose of the matrix to be loaded in Tetrad
##write(t(new.mat1),file="test_conservative_pc_data1.txt",ncolumns=p)

##check the output with Tetrad

amat.tetrad1 <- matrix(0,p,p)
amat.tetrad1[1,6] <- 1
amat.tetrad1[2,4] <- amat.tetrad1[2,6] <- 1
amat.tetrad1[3,4] <- 1
amat.tetrad1[4,6] <- amat.tetrad1[4,9] <- 1

correctEst1 <- all(dag1.amat == amat.tetrad1)
if (!correctEst1) stop("Test sample conservative PC wrong: example 1!")
showProc.time()


## Example 2
p <- 12
n <- 5000
set.seed(p*67+15673)
g <- randomDAG(p,2/(p-1))

## generate n samples of DAG using standard normal error distribution
##random data
load(system.file("external", "test_conservative_pc_data2.rda", package = "pcalg"))
## set.seed(67*67)
## new.mat <- rmvnorm(n,mean=rep(0,p),sigma=trueCov(g))
## save(new.mat, file = "/u/kalisch/research/packages/pcalg/inst/external/test_conservative_pc_data2.rda")
## load(file = "/u/kalisch/research/packages/pcalg/inst/external/test_conservative_pc_data2.rda")

suffStat.data <- list(C=cor(new.mat2),n=n)
## pcAlgo conservative sample
dag2 <- pc(suffStat.data, gaussCItest, alpha=0.005, p=p,
           u2pd="relaxed", conservative=TRUE)

##adjacency matrix
dag2.amat <- as(dag2@graph,"matrix")

##always save the transpose of the matrix to be loaded in Tetrad
##write(t(new.mat),file="test_conservative_pc_data2.txt",ncolumns=p)

##check the output with Tetrad

amat.tetrad2 <- matrix(0,p,p)
amat.tetrad2[1,2] <- amat.tetrad2[1,7] <- 1
amat.tetrad2[2,1] <- amat.tetrad2[2,11] <- 1
amat.tetrad2[3,9] <- amat.tetrad2[3,11] <- amat.tetrad2[3,12] <- 1
amat.tetrad2[4,5] <- amat.tetrad2[4,11] <- amat.tetrad2[4,12] <- 1
amat.tetrad2[5,4] <- amat.tetrad2[5,8] <- amat.tetrad2[5,9] <- 1
amat.tetrad2[7,1] <- 1
amat.tetrad2[8,5] <- amat.tetrad2[8,12] <- 1
amat.tetrad2[10,11] <- 1

correctEst2 <- all(dag2.amat == amat.tetrad2)
if (!correctEst2) stop("Test sample conservative PC wrong: example 2!")
showProc.time()


##PC algorithm population
##_____________________________________________________________________________
## Example 4
p <- 15
set.seed(15673)
g <- randomDAG(p, 2/(p-1))

##population version
suffStat <- list(C=cov2cor(trueCov(g)),n=10^9)
dag4 <- pc(suffStat, gaussCItest, alpha=0.9999, p=p, u2pd="relaxed")
dag5 <- pc(suffStat, gaussCItest, alpha=0.9999, p=p, u2pd="relaxed",conservative=TRUE)

##adjacency matrix
dag4.amat <- as(dag4@graph,"matrix")
dag5.amat <- as(dag5@graph,"matrix")

correctEst4 <- all(dag4.amat == dag5.amat)
if (!correctEst4) stop("Test population conservative PC wrong: example 4!")
showProc.time()


## Example 5
p <- 25
set.seed(1589873)
g <- randomDAG(p,2/(p-1))

##population version
suffStat <- list(C=cov2cor(trueCov(g)),n=10^9)
dag6 <- pc(suffStat, gaussCItest, alpha=0.9999, p=p, u2pd="relaxed")
dag7 <- pc(suffStat, gaussCItest, alpha=0.9999, p=p, u2pd="relaxed",conservative=TRUE)

##adjacency matrix
dag6.amat <- as(dag6@graph,"matrix")
dag7.amat <- as(dag7@graph,"matrix")

correctEst5 <- all(dag6.amat == dag7.amat)
if (!correctEst5) stop("Test population conservative PC wrong: example 5!")
showProc.time()

## Example 6
p <- 35
set.seed(78673)
g <- randomDAG(p,2/(p-1))

##population version
suffStat <- list(C=cov2cor(trueCov(g)),n=10^9)
dag8 <- pc(suffStat, gaussCItest, alpha=0.9999, p=p, u2pd="relaxed")
dag9 <- pc(suffStat, gaussCItest, alpha=0.9999, p=p, u2pd="relaxed", conservative=TRUE)

##adjacency matrix
dag8.amat <- as(dag8@graph,"matrix")
dag9.amat <- as(dag9@graph,"matrix")

correctEst6 <- all(dag8.amat == dag9.amat)
if (!correctEst6) stop("Test population conservative PC wrong: example 6!")

showProc.time()
