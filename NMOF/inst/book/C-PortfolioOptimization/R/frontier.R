# frontier.R -- version 2010-11-24
require(quadprog)
require(NMOF)

# get data
R <- fundData[1L:100L,1L:50L]

# true values and settings
na <- dim(R)[2L]   # number of assets 
ns <- 50L
m  <- colMeans(R)
Sigma <- cov(R)
wsup  <- 1.0       # maximum holding size
winf  <- 0.0       # minimum holding size

# compute frontier
nFP       <- 100L  # number of frontier points
lambdaSeq <- seq(0.001, 0.999, length = nFP) 
A <- array( 1, dim = c(1L,na))
B <- rbind(-diag(na),diag(na))
a <- 1
b <- rbind(array(-wsup, dim = c(na,1L)), 
           array( winf, dim = c(na,1L)))

# matrix for effcient portfolios
pMat <- array(NA, dim = c(na,nFP)) 
for(lambda in lambdaSeq) {
    result <- solve.QP(Dmat = 2*(1-lambda)*Sigma,
                       dvec = lambda*m,
                       Amat = t(rbind(A,B)),
                       bvec = rbind(a,b),
                       meq  = 1L)
    pMat[ ,which(lambda==lambdaSeq)] <- result$solution
}

# plot results
barplot(pMat, legend.text = TRUE, space = 0)