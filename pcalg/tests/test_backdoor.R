library(pcalg)
(doExtras <- pcalg:::doExtras())

## NB: add tests in addition to the simple ones from Maathuis and Colombo (2013)
## in ../man/backdoor.Rd
##    ~~~~~~~~~~~~~~~~~~

`%w/o%` <- function(x, y) x[!x %in% y] #--  x without y
## slightly faster:
`%w/o%` <- function(x, y) x[!match(x, y, nomatch = 0L)]

###-------- DAG ----------------------
set.seed(47)
p <- if(doExtras) 17 else 12
myDAG <- randomDAG(p, prob = 1/4) ## true DAG

## Extract the adjacency matrix of the true DAG
true.amat <- (amat <- as(myDAG, "matrix")) != 0 # TRUE/FALSE <==> 1/0
print.table(1*true.amat, zero.=".") # "visualization"

nodes <- 1:p; names(nodes) <- nodes
cat("Time for many backdoor(.., \"dag\") s : ", system.time(
LL <- lapply(nodes, function(i)
	     lapply(nodes %w/o% i,
		    backdoor,
		    amat = true.amat, x = i, type="dag"))
), "\n")

if(doExtras) {
    for(i in nodes[1:3]) ## Nodes 1,2,3 are all "root" nodes:
	stopifnot(vapply(LL[[i]], identical, NA, y=integer(0)))

    str(LL[-(1:3)]) ## Martin: interesting.. Q: is "this" known? A: yes, basically
} else {
    str(LL)
}

###-------- CPDAG --------------------

## estimate the true CPDAG
myCPDAG <- dag2cpdag(myDAG)
## Extract the adjacency matrix of the true CPDAG
CP.amat <- (as(myCPDAG, "matrix") != 0) # 1/0 <==> TRUE/FALSE

cat("Time for many backdoor(.., \"cpdag\") s : ", system.time(
L2 <- lapply(nodes, function(i)
	     lapply(nodes %w/o% i,
		    backdoor,
		    amat = CP.amat, x = i, type="cpdag"))
), "\n")

str(L2)


###-------- PAG ----------------------

## define nodes 2 and 6 to be latent variables
L <- c(2,6)

## compute the true covariance matrix of g
cov.mat <- trueCov(myDAG)
## transform covariance matrix into a correlation matrix
true.corr <- cov2cor(cov.mat)

## find PAG
## as dependence "oracle", we use the true correlation matrix in
## gaussCItest() with a large "virtual sample size" and a large alpha:
true.pag <- dag2pag(suffStat = list(C = true.corr, n = 10^9),
                    indepTest = gaussCItest,
                    graph=myDAG, L=L, alpha = 0.9999)
PAG.amat <- true.pag@amat # {0 1 2 3}

cat("Time for many backdoor(.., \"pag\") s : ", system.time(
L3 <- lapply(nodes, function(i)
	     lapply(nodes %w/o% i,
		    backdoor,
		    amat = PAG.amat, x = i, type="pag"))
), "\n")

str(L3)


###-------- MAG ----------------------


## find a valid MAG such that no additional edges are directed into
(MAG.amat <- pag2magAM(PAG.amat, 4))

cat("Time for many backdoor(.., \"mag\") s : ", system.time(
L4 <- lapply(nodes, function(i)
	     lapply(nodes %w/o% i,
		    backdoor,
		    amat = MAG.amat, x = i, type="mag"))
), "\n")
## actually this is the *fastest* of the cases

str(L4)
