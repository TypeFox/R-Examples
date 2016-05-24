factorQR <- function(factorForm, nonFactorForm=NULL, dataSet=parent.frame(), pQuant = 0.5, whichFactor = NULL, 
nSamp=5000, burn=0, thin = 0, cTau0=1, dTau0=1, cPsi0=1, dPsi0=1,
sig0=1, mu0=1, R0=NULL, nu0=NULL, B0s=0.01, B0Beta = .01,
betaZero = NULL, 
PhiZero = NULL, invPsiZero = NULL, LambdaZero=NULL, LambdaSZero=NULL, OmegaZero = NULL, verbose=TRUE, storeOmega=FALSE, 
latentInteract=FALSE, interactX=NULL, whichFactorInteract=NULL){

whichFactorInteract <- whichFactorInteract - 1
	
## need to deal with two formulas: 
captureFormula <- function(formula, data=dataSet){
        call <- match.call()
        mf <- match.call(expand = FALSE)
        mf[[1]] <- as.name("model.frame.default")
        mf <- eval.parent(mf)
        return(mf)
    }

interManInd <- !is.null(interactX)
if(latentInteract && interManInd)
stop("This function does not currently support latent/latent interactions at the same time as latent/manifest interations.")

numManInd <- length(whichFactorInteract)

mf1 <- captureFormula(factorForm, data=dataSet)
Terms <- attr(mf1, "terms")
y=model.response(mf1)
x=model.matrix.default(Terms, mf1)
x <- x[,apply(x,2,sd) > 0] # remove constant column, if it exists
x <- scale(x, scale=FALSE)  #center x
xnames <- colnames(x)


if(!is.null(nonFactorForm)){
mf2 <- captureFormula(nonFactorForm, data=dataSet)
Terms <- attr(mf2, "terms")
nonFactorX <- model.matrix.default(Terms, mf2)
## center X: 
#nonFactorX[,apply(nonFactorX, 2, sd) > 0] <- scale(nonFactorX[, apply(nonFactorX, 2, sd)>0], scale=FALSE)
betaNames <- colnames(nonFactorX)
}
else {nonFactorX <- NULL
	betaNames <- NULL}


call <- match.call()

xLength <- length(x[1,])
xBetLen <- length(nonFactorX[1,])


if(is.null(whichFactor))
whichFactor <- rep(1, xLength)

if(length(whichFactor) != xLength)
stop("`whichFactor' must have length equal to the number of columns in `x' matrix")

whichFactor <- as.numeric(as.factor(whichFactor))

# function that keeps track of dimensions of Lambda:
indexFcn <- function(x){
	return(whichFactor[x])
	}	
		
## Find variables that have two levels: 		
whichDichot <- as.numeric(unlist(lapply(apply(x, 2, unique), length)) == 2)
		
#if(is.null(whichDichot))
#whichDichot <- rep(0, xLength)

#if(length(whichDichot) != xLength)
#stop("`whichDichot' must have length equal to the number of columns in `x' matrix")

## add indicator of number of dichotomous variables
nDichot <- sum(whichDichot > 0)
dichotFcn <- (1:xLength)[whichDichot>0]
#whichDichot <- c(whichDichot, nDichot)



		
		
nObs <- dim(x)[1]
nFact <- length(levels(as.factor(whichFactor)))
if(is.null(PhiZero)) PhiZero <- diag(rep(1, nFact))
Phi <- PhiZero
invPhi <- solve(Phi)
		
if(!is.null(interactX)){
	if(is.null(dim(interactX))){
		stop("`interactX' must be a matrix or data frame.")}
		if(dim(interactX)[1] != nObs){
		stop("`interactX' must be a matrix or data frame with the number of rows equal to the number of observations.")}
		if(dim(interactX)[2] != numManInd){
		stop("`interactManifest' must be a matrix or data frame with the number of columns equal to the length of `whichFactorInteract'. ")}
		}		
		
  if ((burn < 0) || (burn %% 1 !=0))
    stop("`burn' should be a non-negative integer.") 
  if (thin < 0 || (thin %% 1 !=0))
    stop("`thin' should be a non-negative integer.")

  keep <- thin + 1
		


if(is.null(invPsiZero)) invPsiZero <- 1
if(length(invPsiZero) == 1) invPsiZero <- rep(invPsiZero, xLength)
if(length(invPsiZero) != xLength) stop("The length of `invPsiZero' should be 1 or the number of columns in `x'")

invPsi <- invPsiZero



if(is.null(OmegaZero))
Omega <- matrix(rnorm(nObs*nFact), nr=nObs)
else {
	if(length(OmegaZero) != nObs*nFact){
		cat("Warning: OmegaZero is of the wrong length.  Generating new OmegaZero.\n")
		Omega <- matrix(rnorm(nObs*nFact), nr=nObs)
		}
	else
	Omega <- OmegaZero
}


if(!latentInteract && !interManInd){
if(!(length(B0s) %in% c(1,nFact, nFact^2))) stop("`B0s' must be a vector or square matrix with dimension equal to the number of factors plus one or a scalar.")
if(length(B0s) == 1)
B0s <- diag(rep(B0s, nFact), nrow=nFact)
else if(length(B0s) == nFact)
B0s <- diag(B0s, nrow=nFact)
}

if(latentInteract){
if(!(length(B0s) %in% c(1,nFact+1, (nFact+1)^2))) stop("`B0s' must be a vector or square matrix with dimension equal to the number of factors or a scalar.")
if(length(B0s) == 1){
B0s <- diag(rep(B0s, nFact + 1), nrow=(nFact+1))
}
else if(length(B0s) == (nFact+1))
B0s <- diag(B0s, nrow=(nFact+1))
}

if(interManInd){
if(!(length(B0s) %in% c(1,nFact+numManInd, (nFact+numManInd)^2))) stop("`B0s' must be a vector or square matrix with dimension equal to the number of factors plus factor interaction terms or a scalar.")
if(length(B0s) == 1){
B0s <- diag(rep(B0s, nFact + numManInd), nrow=(nFact+numManInd))
}
else if(length(B0s) == (nFact+numManInd))
B0s <- diag(B0s, nrow=(nFact+numManInd))
}






if(!is.null(nonFactorX)){

if(!(length(B0Beta) %in% c(1,xBetLen, xBetLen^2))) stop("`B0Beta' must be a vector or square matrix with dimension equal to the number of factors or a scalar.")

if(length(B0Beta) == 1)
B0Beta <- diag(rep(B0Beta, xBetLen))

else if(length(B0Beta) == xBetLen)
B0Beta <- diag(B0Beta)
}


if(is.null(nonFactorX)){
	B0Hold <- B0s
	}
else{
	if(!latentInteract){
	B0Hold <- matrix(0, nr=nFact + xBetLen + numManInd, nc=nFact + xBetLen + numManInd)
	B0Hold[1:(nFact+numManInd), 1:(nFact+ numManInd)] <- B0s
	B0Hold[(nFact+numManInd + 1):(nFact+numManInd+xBetLen), (nFact+numManInd + 1):(nFact+numManInd+xBetLen)] <- B0Beta
	}
	
	if(latentInteract){
	B0Hold <- matrix(0, nr=nFact + xBetLen + 1, nc=nFact + xBetLen + 1)
	B0Hold[1:(nFact+1), 1:(nFact + 1)] <- B0s
	B0Hold[(nFact + 2):(nFact+xBetLen +1), (nFact + 2):(nFact+xBetLen + 1)] <- B0Beta
		}
	}

	
if(is.null(R0)){
	R0 <- diag(rep(1,nFact))
	}	
else if(length(R0) == 1) R0 <- diag(rep(R0, nFact))
else if(length(R0) == nFact) R0 <- diag(R0)
else if(length(R0) != nFact^2) 	stop("`R0' must be a vector or square matrix with dimension equal to the number of factors or a scalar.")
	

if(is.null(nu0)) nu0 <- nFact + 1

if(nu0 < (nFact)){
	cat("Warning: nu0 changed to ", nFact, " for a proper prior.\n", sep="")
	nu0 <- nFact
	}		



if(!is.null(betaZero))
beta <- betaZero
else beta <- rep(0, xBetLen)

if(is.null(LambdaZero))
LambdaZero <- 1

if(length(LambdaZero) == 1)
LambdaZero <- rep(LambdaZero, xLength)

if(length(LambdaZero) != xLength)
stop("LambdaZero should be NULL, a scalar, or have length equal to the number of columns in `x'.")


Lambda <- matrix(0, nr=xLength, nc=nFact)
for(i in 1:xLength){
Lambda[i,indexFcn(i)] <- 1
}

fixedLambdaRows1 <- (1:xLength)[!duplicated(Lambda)]
fixedLambdaRows <- as.numeric(!duplicated(Lambda))

for(i in 1:xLength){
	if(!(i %in% fixedLambdaRows1))
	Lambda[i,indexFcn(i)] <- LambdaZero[i]
	}

if(is.null(LambdaSZero))
LambdaSZero <- 0.1

if(length(LambdaSZero) == 1)
LambdaSZero <- rep(LambdaSZero, nFact)

### if there are zeros in the response, starting LambdaSZero =0 will crash it. ###
if(sum(LambdaSZero == 0) == length(LambdaSZero)){
	LambdaSZero <- LambdaSZero + rnorm(length(LambdaSZero), mean=0, sd=.01)
	cat("Warning: Starting with LambdaSZero at zero can cause problems in the sampler.  Noise with sd=0.01 added to LambdaSZero.")
	}

if(length(LambdaSZero) != nFact)
stop("`LambdaSZero' should be NULL, a scalar, or a vector with length equal to the number of factors in the model.")	
	
Lambda <- rbind(Lambda, LambdaSZero)

if(sum((LambdaZero[fixedLambdaRows1] - 1)^2) != 0)
cat("Warning.  Components of `LambdaZero' corresponding to founding factors have been changed to 1.")


  if (verbose) {
  	cat("There are ", nFact, " factors.\n", sep="")
    cat("The dimension of beta is ", xBetLen, ".\n", sep="")
    cat("The are ", nObs, " observations.\n\n", sep="")
  } 



### start sampler ###
# n.par <- (Lambda_{-q} + Lambda_q + Psi + Phi + beta + tau + mu)
n.par <- 2*xLength + nFact + nFact * (nFact + 1)/2 +  xBetLen + 1 + nDichot

if(latentInteract) n.par <- n.par + 1

if(interManInd) n.par <- n.par + numManInd

if(storeOmega) n.par <- n.par + nObs * nFact

	
indexFcn2 <- indexFcn(1:xLength)	

## need place holders even if nonFactorX = NULL:
if(xBetLen == 0){
	
	nonFactorX=0
	beta=0
	xBetLen=0
	}
	
	
if(is.null(interactX))
interactX <- rep(0, nObs)	

	lamNames <- NULL
	for(i in 1:nFact) lamNames <- c(lamNames, paste("RespLambda-", levels(as.factor(whichFactor))[i], sep=""))
	if(latentInteract) lamNames <- c(lamNames,"FactorInteract")
	if(interManInd){
		for(i in 1:numManInd){
			lamNames <- c(lamNames, paste(lamNames[whichFactorInteract[i+1]], "-Interact-", i, sep=""))
			}
		}
	for(i in 1:xLength) lamNames <- c(lamNames, paste("Lambda_{",xnames[i], ",", indexFcn(i), "}", sep=""))

	
	phiNames <- NULL
	for(j in 1:nFact)
	for(k in 1:nFact)
	if(j<=k) phiNames <- c(phiNames, paste("Phi_{", j, ",", k,"}", sep=""))
	
	psiNames <- NULL
	for(j in 1:xLength) psiNames <- c(psiNames, paste("Psi_", j, sep=""))
	
	muNames <- NULL
	if(nDichot > 0){
	for(i in 1:nDichot){
		muNames <- c(muNames, paste("mu_{", xnames[dichotFcn[i]],"}", sep=""))
	}}
	
#	for(i in 1:xBetLen) betaNames <- c(betaNames, paste("beta", i, sep=""))
	
	nameVec <- c(lamNames, betaNames, phiNames, psiNames, "tau", muNames)

if(verbose)
	cat("Starting sampler... \n")
	
#	return(list(nObs=nObs, nSamp=nSamp, xLength=xLength, nFact=nFact, xBetLen=xBetLen, x=x, y=y, nonFactorX=nonFactorX, pQuant=pQuant, which.factor=which.factor, indexFcn2=indexFcn2, fixedLambdaRows=fixedLambdaRows, burn = burn, thin=thin, cTau0=cTau0, dTau0=dTau0, betaPsiZero0=betaPsiZero, alphaPsiZero=alphaPsiZero, sig0=sig0, mu0=mu0, R0=R0, nu0=nu0, B0Hold=B0Hold, beta=beta, invPhi=invPhi, invPsi=invPsi, Lambda=Lambda, Omega=Omega))

#if(nDichot == 0){

param <- .C("cQuantRegFact",as.integer(nObs), as.integer(nSamp), as.integer(xLength), as.integer(nFact), 
	as.integer(xBetLen), as.double(x), as.double(y), as.double(nonFactorX), 
	as.double(pQuant), as.integer(nDichot), 
	as.integer(whichDichot), as.integer(dichotFcn), as.integer(indexFcn2), as.integer(fixedLambdaRows), 
	as.integer(burn), as.integer(thin),
	as.double(cTau0), as.double(dTau0), as.double(dPsi0), 
	as.double(cPsi0), as.double(sig0), as.double(mu0),
	as.double(R0), as.integer(nu0), as.double(B0Hold), as.double(beta), 
	as.double(invPhi), as.double(invPsi),
	as.double(Lambda), as.double(Omega), as.integer(verbose), 
	as.integer(storeOmega),
	as.integer(latentInteract), as.double(interactX), as.integer(numManInd), as.integer(whichFactorInteract),
	pdStore = double(n.par*floor((nSamp-burn)/keep)), 
	PACKAGE="factorQR")$pdStore
	
#	}
	

param <- matrix(param, ncol = n.par, nrow = floor((nSamp-burn)/keep), byrow=TRUE)

if(storeOmega){
	omega <- param[,-(1:(n.par-nObs * nFact))]
	param <- param[,1:(n.par-nObs * nFact)]
	}
else omega <- NULL	

dimnames(param)[[2]] <- nameVec

nReg = xBetLen + nFact + xLength + numManInd
if(latentInteract) nReg <- nReg + 1

nReg <- nReg + numManInd

nFactInt <- nFact
if(latentInteract) nFactInt <- nFactInt + 1

nFactInt <- nFactInt + numManInd

res <- list(param = param, call = call, nReg = xBetLen + nFact + xLength, betLen = xBetLen, nObs=nObs, burn = burn, thin=thin, nSamp=nSamp, nFact=nFact, nFactorX = xLength, omega=omega, nFactInt=nFactInt)

class(res) <- "factorQR"

return(res)
}


