makeData <- function(N, whichFactor, pQuant=.5, Lambda = NULL, LambdaQ = NULL, Phi = NULL, lapScale=1, Psi = NULL, interact=0){


## N: sample size
## pQuant: pth quantile
## Lamb: Lambda_{-q}
## LambQ: Lambda_q
## Phi: Phi (covariance structure for latent factors)
## lapScale: scale of Laplace error
## Psi: Psi (error variances for epsilon_{-q})

## ralap taken from VGAM package under GPL-2
## See cran.r-project.org/web/packages/VGAM/index.html
ralap <- function (n, location = 0, scale = 1, tau = 0.5, kappa = sqrt(tau/(1 - 
    tau))) 
{
	
is.Numeric <-  function (x, allowable.length = Inf, integer.valued = FALSE, 
    positive = FALSE) 
if (all(is.numeric(x)) && all(is.finite(x)) && (if (is.finite(allowable.length)) length(x) == allowable.length else TRUE) && (if (integer.valued) all(x == 
    round(x)) else TRUE) && (if (positive) all(x > 0) else TRUE)) TRUE else FALSE
	
    use.n = if ((length.n <- length(n)) > 1) 
        length.n
    else if (!is.Numeric(n, integ = TRUE, allow = 1, posit = TRUE)) 
        stop("bad input for argument 'n'")
    else n
    location = rep(location, len = use.n)
    scale = rep(scale, len = use.n)
    tau = rep(tau, len = use.n)
    kappa = rep(kappa, len = use.n)
    ans = location + scale * log(runif(use.n)^kappa/runif(use.n)^(1/kappa))/sqrt(2)
    indexTF = (scale > 0) & (tau > 0) & (tau < 1) & (kappa > 
        0)
    ans[!indexTF] = NaN
    ans
}



xLength <- length(whichFactor)

indexFcn <- function(x){
	return(whichFactor[x])
	}	
	
nFact <- length(levels(as.factor(whichFactor)))

if(is.null(Lambda)) Lambda <- rep(1, xLength)
if(is.null(LambdaQ)) LambdaQ <- rep(0, nFact)
if(is.null(Phi)) Phi <- diag(rep(1, nFact))
if(is.null(Psi)) Psi <- diag(rep(1, xLength))

lambda <- matrix(0, nr=xLength, nc=nFact)
for(i in 1:xLength){
lambda[i,indexFcn(i)] <- 1
}

fixedLambdaRows1 <- (1:xLength)[!duplicated(lambda)]
fixedLambdaRows <- as.numeric(!duplicated(lambda))

for(i in 1:xLength){
	if(!(i %in% fixedLambdaRows1))
	lambda[i,indexFcn(i)] <- Lambda[i]
	}
lambda <- rbind(lambda, LambdaQ)

cat("Lambda is: \n\n")
print(lambda)

rtPhi <- chol(Phi)

Omega <- matrix(rnorm(N*nFact), nc=nFact) %*% rtPhi

xMat <- t(lambda %*% t(Omega))
psi <- diag(Psi)
invPsi <- 1/psi
trueInvPsi <- invPsi
rtPsi <- chol(Psi)

if(interact != 0)
yInteract <- Omega[,1] * Omega[,2] * interact
else yInteract <- 0

errorMat <- matrix(rnorm(N*xLength), nc=xLength) %*% rtPsi

yError <- ralap(N, location = 0, scale = lapScale, tau=pQuant)
yVec <- xMat[,xLength+1]
yVec <- yVec + yError + yInteract


xMat <- xMat[,-(xLength + 1)]
xMat <- xMat + errorMat
dataSet <- data.frame(yVec, xMat)
colNames <- "Y"
for(i in 1:xLength) colNames <- c(colNames, paste("X", i, sep=""))
names(dataSet) <- colNames
return(dataSet)
}
