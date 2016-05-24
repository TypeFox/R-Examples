bayesQR <- function(formula, dataSet=parent.frame(), pQuant = 0.5,  
nSamp=5000, burn=0, thin = 0, C0=1, D0=1, 
B0=0.01, betaZero = NULL,  verbose=TRUE){
## x: matrix of modeled x components
## y: vector of responses
## nonFactorX: matrix of standard covariates -- not in RHS -- 
##     to model y; **does not include intercept--add if needed**
## pQuant: quantile to model
## which.factor: vector of indicators to show factor grouping
## 		e.g., if which.factor = c(1,1,1,2,2,2), it would mean the first three 
## 		dimensions of x group together and the last three group as well.
## nSamp: number of MCMC iterations
## burn: iterations of burn-in
## thin: number of iterations to skip between stored values
## c0: hyperparameter for tau (shape)
## d0: hyperparameter for tau (scale)
## betaPsiZero: hyperparameter for prior for Psi (rate)
## alphaPsiZero: hyperparameter for prior for Psi (shape)
## sig0: hyperparameter for Lambda parameters: prior variance is psi_j*sig0 
## mu0: prior mean for Lambda_{-q} components
## R0: hyperparameter for Phi's prior -- actually R_0^{-1} from the paper; Identity matrix default
## nu0: prior dof for inverse Wishart
## B0q: prior for Lambda_q: prior precision (centered at zero)
## which.categorical: matrix of indicators for which x variables are indicators (coded 0, 1)
## store: if QOnly, only stores qth row of Lambda; if Psi=TRUE, stores Psi values; etc.
## betaZero: initial value for beta (corresponding to nonFactorX)
## PhiZero: initial value of Phi
## invPsiZero: initial value of Psi^{-1}
## lambdaZero: inital value of Lambda
## printProgress:  if not null, prints update every printProgress iterations
	

#storeCount <- 1

captureFormula <- function(formula, data=dataSet){
        call <- match.call()
        mf <- match.call(expand = FALSE)
        mf[[1]] <- as.name("model.frame.default")
        mf <- eval.parent(mf)
        return(mf)
    }


mf1 <- captureFormula(formula, data=dataSet)
Terms <- attr(mf1, "terms")
y=model.response(mf1)
x=model.matrix.default(Terms, mf1)
betaNames <- colnames(x)


call <- match.call()

xLength <- length(x[1,])
		
		
nObs <- dim(x)[1]		
		
  if ((burn < 0) || (burn %% 1 !=0))
    stop("`burn' should be a non-negative integer.") 
  if (thin < 0 || (thin %% 1 !=0))
    stop("`thin' should be a non-negative integer.")

  keep <- thin + 1
		



if(!(length(B0) %in% c(1,xLength, xLength^2))) stop("`B0' must be a vector or square matrix with dimension equal to the number of columns in the design matrix or a scalar.")

if(length(B0) == 1){
B0 <- diag(rep(B0, xLength), nrow=xLength)  #need nr=xLength here... otherwise messes up with nr=1
#print(B0q)
}




if(!is.null(betaZero))
	beta <- betaZero
else beta <- 0
	
if(length(beta) == 1)
beta <- rep(beta, xLength)	
	

#print(dim(x))
#cat("dim of x\n")
#print(beta)

XBeta <- x %*% matrix(beta, nc=1)

if(sum(XBeta == y) > 0){
	cat("Warning.  Some starting values fit the data perfectly, which is not allowed.  Some noise has been added.\n")
	}
	
XBeta[XBeta == y] <- XBeta + sample(c(-.01, .01), 1)




  if (verbose) {
    cat("The dimension of beta is ", xLength, ".\n", sep="")
    cat("The are ", nObs, " observations.\n\n", sep="")
  } 


n.par <- xLength + 1


	


nameVec <- c(betaNames, "tau")

if(verbose)
	cat("Starting sampler... \n")
	
param <- .C("cBayesQuantReg",as.integer(nObs), as.integer(nSamp), as.integer(xLength),  
	as.double(x), as.double(y), as.double(XBeta), 
	as.double(pQuant), 
	as.integer(burn), as.integer(thin),
	as.double(C0), as.double(D0), 
	as.double(B0), as.double(beta), 
	as.integer(verbose), 
	pdStore = double(n.par*floor((nSamp-burn)/keep)), 
	PACKAGE="factorQR")$pdStore

param <- matrix(param, ncol = n.par, nrow = floor((nSamp-burn)/keep), byrow=TRUE)

dimnames(param)[[2]] <- nameVec

res <- list(param = param, call = call, betLen = xLength, nObs=nObs, burn = burn, thin=thin, nSamp=nSamp)

class(res) <- "bayesQR"

return(res)
}


