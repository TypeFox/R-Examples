library(pcalg)
source(system.file(package="Matrix", "test-tools-1.R", mustWork=TRUE))
##--> showProc.time(), assertError(), relErrV(), ...
R.home(); sessionInfo() # helping package maintainers to debug ...
.libPaths()
packageDescription("pcalg")
packageDescription("Matrix")

## load the functions for the simulations of this paper
## source("/u/colombo/Diss/RAusw/First_paper_RFCI/functions_for_the_simulations.R")

## RFCI improves the output
##______________________________________________
## Input:  L1=1; L2=2; X1=6; X2=4; X3=3; X4=5; X5=7; X6=8
## Output: X1=4; X2=2; X3=1; X4=3; X5=5; X6=6

amat <- t(matrix(c(0,0,0,0,0,0,0,0, 0,0,0,0,1,1,0,0,
                   0,1,0,1,0,1,0,0, 1,0,0,0,0,1,0,0,
                   0,0,0,0,0,1,0,0, 1,1,0,0,0,0,0,0,
                   0,0,0,1,1,0,0,0), 8,8))
colnames(amat) <- rownames(amat) <- as.character(1:8)
Matrix::Matrix(amat) # to "visualize"
L <- c(1,2)
V <- as.character(1:8)
edL <- setNames(vector("list", length=length(V)), V)
edL[[6]] <- list(edges=NULL, weights=NULL)
edL[[8]] <- list(edges=NULL, weights=NULL)
edL[[4]] <- list(edges=c(7,8),   weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[3]] <- list(edges=c(4,5,8), weights=c(abs(rnorm(1)),abs(rnorm(1)),abs(rnorm(1))))
edL[[5]] <- list(edges=c(6,8),   weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[7]] <- list(edges= 8,       weights=abs(rnorm(1)))
edL[[1]] <- list(edges=c(4,6),   weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[2]] <- list(edges=c(5,7),   weights=c(abs(rnorm(1)),abs(rnorm(1))))
g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")
if(dev.interactive())
    plot(g)

## Compute the true covariance matrix of g
cov.mat <- trueCov(g)
## Delete rows and columns which belong to L
true.cov <- cov.mat[-L,-L]
## Transform it in a correlation matrix
true.corr <- cov2cor(true.cov)

suffStat <- list(C=true.corr, n=10^9)
showSys.time(pop.fci1 <-
             fci(suffStat, gaussCItest, labels=V[-L],
                 alpha=0.9999, doPdsep=TRUE,verbose=FALSE)@amat)

showSys.time(pop.rfci1 <-
             rfci(suffStat, gaussCItest, labels=V[-L],
                  alpha=0.9999, verbose=FALSE)@amat)

if (any(pop.fci1 != pop.rfci1)) {
  stop("Test of RFCI wrong: small example!")
}

## Thomas' example (version number 8) about discriminating path orientation rule

V <- as.character(1:25)
edL <- setNames(vector("list", length=length(V)), V)
edL[[ 1]] <- list(edges=c(14,18),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[ 2]] <- list(edges=c(16,18),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[ 3]] <- list(edges=c(16,24),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[ 4]] <- list(edges=c(18,24),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[ 5]] <- list(edges=c(15,25),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[ 6]] <- list(edges=c(17,19),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[ 7]] <- list(edges=c(14,19),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[ 8]] <- list(edges=c(14,24),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[ 9]] <- list(edges=c(19,20),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[10]] <- list(edges=c(20,25),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[11]] <- list(edges=c(23,25),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[12]] <- list(edges=c(22,24),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[13]] <- list(edges=c(21,23),weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[14]] <- list(edges=NULL,    weights=NULL)
edL[[15]] <- list(edges=c(16,17,24),weights=c(abs(rnorm(1)),abs(rnorm(1)),abs(rnorm(1))))
edL[[16]] <- list(edges=c(19,25),   weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[17]] <- list(edges=c(18,24,25),weights=c(abs(rnorm(1)),abs(rnorm(1)),abs(rnorm(1))))
edL[[18]] <- list(edges=c(21,25),   weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[19]] <- list(edges=c(23,24,25),weights=c(abs(rnorm(1)),abs(rnorm(1)),abs(rnorm(1))))
edL[[20]] <- list(edges= 24,        weights=abs(rnorm(1)))
edL[[21]] <- list(edges=c(22,25),   weights=c(abs(rnorm(1)),abs(rnorm(1))))
edL[[22]] <- list(edges= 25,        weights=abs(rnorm(1)))
edL[[23]] <- list(edges= 24,        weights=abs(rnorm(1)))
edL[[24]] <- list(edges=NULL,weights=NULL)
edL[[25]] <- list(edges=NULL,weights=NULL)
(g <- new("graphNEL", nodes=V, edgeL=edL,edgemode="directed"))

if(dev.interactive())
    plot(g)

## Latent variables (all having no parents):
L <- c(1:13)

## Compute the true covariance matrix of g
cov.mat <- trueCov(g)
## Delete rows and columns which belong to L
true.cov <- cov.mat[-L,-L]
## Transform it in a correlation matrix
true.corr <- cov2cor(true.cov)
suffStat <- list(C=true.corr, n=10^9)
p.tr <- dim(true.corr)[1]
showSys.time(pop.fci2  <-  fci(suffStat, gaussCItest, p=p.tr,
                               alpha=0.9999, doPdsep=TRUE)@amat)
showSys.time(pop.rfci2 <- rfci(suffStat, gaussCItest, p=p.tr,
                               alpha=0.9999)@amat)

if (any(pop.fci2 != pop.rfci2)) {
  stop("Test of RFCI wrong: big example!")
}
