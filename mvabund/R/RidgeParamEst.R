################################################################################
# RIDGEPARAMEST Maximum likelihood estimation of the ridge parameter
# by cross-validation
################################################################################

ridgeParamEst<-function(dat, X, weights = rep(1,times=nRows), refs, 
	tol=1.0e-010, only.ridge=FALSE,  doPlot=FALSE, col="blue", type="l", ...) {
	
 dat <- as.matrix(dat)
 X   <- as.matrix(X)
 nRows   <- nrow(dat)
 if (nRows==1 | nRows != nrow(X)) {
    stop("The row-number of your data is not correct.")
 }

 nVars     <- ncol(dat)
 nCols		 <- ncol(X)
 miss.refs <- missing(refs)
 if (miss.refs) { 
	# Set REFS, which determines the method of cross-validation.
 
    if (nRows<=20)    {
        nFolds      <- nRows
    } else if (nRows<=40) {
        nFolds      <- 10
    } else    {
      
	  nFolds      <- 5
    }
	
    nGroupI         <- ceiling(nRows/nFolds) # max number in each CV group.
    refs            <- rep.int( 1:nFolds,times= nGroupI )
    # length>nRows, but next line also ensures that the length is the same again.
    # Put refs into a random order.
    refs            <- refs[ order( runif(nRows) )  ]
 } else {
  # Set default values for other input arguments.
	if(length(refs)!=nRows)
    stop("the length of 'refs' must be the sample size of 'dat'")
    nFolds <- length( unique(refs) )
}
isSing     <- FALSE
mnPred     <- matrix(0,nrow=nRows,ncol=nVars)

sdInvTrain <- RTrain <- list()
sdInvTr    <- eyeTr  <- matrix(0, nrow=nVars, ncol=nVars)
eyeTr[1+0:(nVars-1)*(nVars+1)] <- 1 

# loop
singular <- 0
for (iFold in 1:nFolds) { 
    isTrain   <- (refs != iFold)     # index for Training data
	  isTest    <- (refs == iFold)     # index for Test data
    XTrain    <- X[isTrain,, drop = FALSE]    # the training data
   	datTrain  <- dat[isTrain,, drop = FALSE]  # the training design matrix data set
	  qrTrain	  <- qr(XTrain*sqrt(weights[isTrain]), tol=tol)
	  rank	  <- qrTrain$rank
	
	if(qrTrain$rank < nCols & miss.refs & nFolds < nRows){
		return( ridgeParamEst(dat=dat, X=X, weights=weights, refs=1:nRows, 
			tol=tol, only.ridge= only.ridge, doPlot=doPlot, col=col,type=type,...))
	}
	XTrain	<- XTrain[,qrTrain$pivot[1:rank], drop = FALSE ]
	# Obtain the training estimate.
	betaTrain <- chol2inv(qrTrain$qr[1:rank,1:rank, drop=FALSE]) %*% 
		t(XTrain*weights[isTrain]) %*% datTrain
  # Obtain the predicted data with the training beta.
	mnPred[isTest,] <- X[isTest,qrTrain$pivot[1:rank], drop = FALSE] %*% betaTrain 	

	# Obtain the training residuals
  resTrain  <- datTrain - XTrain %*% betaTrain
  s         <- cov(resTrain)   # cov of training residuals

	vr        <- c(s)[1+0:(nVars-1)*(nVars+1)]
  # Variables that should basically be remove from furhther analysis
  isVr0     <- (vr < 10^-10)
	seql.isVr0	 <- 1:sum(isVr0)
	#  vr(isVr0) <- 10^-10 # basically removes these variables from further analyses...
	sdInvTr[1+0:(nVars-1)*(nVars+1)]  <- 1/ sqrt(vr)
	sdInvTrain[[iFold]] <-  sdInvTr
  RTrain[[iFold]] <- sdInvTrain[[iFold]] * s * sdInvTrain[[iFold]]
	# these get var=0, cov=eye
  # Inf values are set to 1
  sdInvTrain[[iFold]][isVr0,isVr0] <- RTrain[[iFold]][isVr0,isVr0] <-
      eyeTr[seql.isVr0,seql.isVr0]  # = diag(sum(isVr0))

	
	rankRT <- qr(RTrain[[iFold]])$rank

  if ( rankRT< nVars) isSing <- TRUE  # flag for k choice.
}

resPred <- dat - mnPred # predicted residuals

# BEGIN iterative estimation of the CV penalized likelihood
    ridgeParameter <- LLSolve(RTrain=RTrain,resPred=resPred,sdInvTrain=sdInvTrain,
		refs=refs,nFolds=nFolds,tol=tol, max.iter=5)
	  
# END iterative estimation

if (!only.ridge) {
   # estimation of the minimal CV penalized likelihood
   # minLL may be complex, if complexDisc set FALSE (default: TRUE)
   minLL <- LLCalc(lambdas=ridgeParameter,isSing=isSing,nFolds=nFolds,resPred=resPred,
   refs=refs,RTrain=RTrain,nVars=nVars,sdInvTrain=sdInvTrain, tol=tol) 

} else {
   minLL<-NULL
}

if (doPlot) {

    lambdaValues <- seq(from=0.01, to =1, by=0.005)
    if (isSing)  {
      lambdaValues<-lambdaValues(lambdaValues!=1) # remove k=1 for singular cases
    }
	
   LLPred <- LLCalc(lambdas=lambdaValues,isSing=isSing,nFolds=nFolds,resPred=resPred,
	  refs=refs,RTrain=RTrain,nVars=nVars,sdInvTrain=sdInvTrain, tol=tol)

   makePlot(lambdas=lambdaValues, LLPred=LLPred, ridgeParameter=ridgeParameter,
    col=col, type=type, ...)
}

  list(ridgeParameter=ridgeParameter, minLL=minLL)

}


################################################################################
# makePlot: function for producing a plot of lambdas and LLPred                #
################################################################################
makePlot <- function (lambdas, LLPred,ridgeParameter, col="blue", type="l",
	xlab=NULL,ylab=NULL,main=NULL,...) {

if (missing(ylab))	ylab  <- "-logL/2"
if (missing(xlab))	xlab  <- expression(lambdas)

plot(lambdas, -LLPred/2, col=col, type=type, ylab=ylab , xlab=xlab,... )
if(missing(main)) {

	if (! missing(ridgeParameter))	{
	
		t <- eval(substitute( expression(hat(lambdas)==u),list(u=ridgeParameter )))
	} else  t <- NULL

	title(t)
}

}


################################################################################
# LLSolve: Using the derivative of the likelihood to estimate ridgeParameter   #
################################################################################
LLSolve <- function(RTrain,resPred,sdInvTrain,refs,nFolds,lambdaInit= 0.5,
iter=1,ridgeParameter,LL=-Inf, tol=1.0e-050, max.iter = 5) {

nVars <- ncol(resPred)
nTest <- matrix( NaN,nrow=nVars,ncol=nFolds)
ZZD   <- LamD <- matrix(nrow=nVars, ncol= nFolds)
for (iFold in 1:nFolds)
  {
     eigRT         <- eigen( RTrain[[iFold]] )
     P             <- eigRT$vectors
     isTest        <- (refs == iFold)
     nTest[,iFold] <- sum(isTest)
     Zi            <- resPred[isTest, ,drop = FALSE] %*% sdInvTrain[[iFold]] %*% P
	   ZZD[,iFold]   <- matrix(1,1,ncol=nrow(Zi)) %*% Zi^2
	   # Obtain eigenvalues of Trainings Corr diag(Lambda)
     LamD[,iFold]  <- eigRT$values
  }
LamD <- Re(LamD)

i        <- 0
lamChange  <- 1
ridgeParameter.new     <- lambdaInit
kappa        <- 1/ridgeParameter.new - 1

while (lamChange>10^(-6) & i<20)
  {
     w.k      <- nTest * ( LamD + kappa )^ (-2)
     num      <- (ZZD/nTest - LamD) * w.k    # componentwise
     kappa      <- sum(num) / sum(w.k)
	   kappa      <- max(kappa, 0) # Ensure that 0<=lambda<=1
     kappa      <- min(kappa, 1/tol) # Ensure that  sum(num) / sum(w.k) is not NaN

     lamChange  <- abs( ridgeParameter.new - ( kappa + 1 ) ^ (-1) )
     ridgeParameter.new     <- ( kappa + 1 ) ^ (-1)
     i          <- i+1

     if (i==20)
     {
         warning("Has not converged in 20 iterations.")
     }
  }

deriv2 <- -1/2 + ZZD / (LamD + kappa)

deriv2   <- sum( deriv2 * w.k)
 
# Set the ridgeParameter.
if (missing(ridgeParameter))
    ridgeParameter <- ridgeParameter.new

if (deriv2<0)
{
  warning("Initial solution is not a maximum. Looking for alternatives...")
    LL_new = LLCalc(lambdas=ridgeParameter,isSing=FALSE,nFolds=nFolds,
	resPred=resPred,refs=refs,RTrain=RTrain,nVars=ncol(resPred),
    sdInvTrain=sdInvTrain, tol=tol)
    if (LL_new < LL)
    {
        ridgeParameter = ridgeParameter.new
        LL = LL_new
    }
    lambdaInit <- lambdaInit*0.5
    iter <- iter+1  # iteration for the next try of LLSolve with changed input

    # max.iter th iteration: stop and present plot with results so far
    if (iter >= max.iter)
    {
        lambdas <- seq(from=0.01, to=0.99, by=0.01 )
        LLPred  <- LLCalc(lambdas=lambdas,isSing=FALSE,nFolds=nFolds,
	           resPred=resPred,refs=refs,RTrain=RTrain,nVars=ncol(resPred),
              sdInvTrain=sdInvTrain, tol=tol)
        makePlot( lambdas, LLPred, ridgeParameter, c(1, 1, 1) )

        if (min(LLPred)<LL)
            stop("No maximum likelihood could be found(!?). See plot.")
        else
            warning("Check plot to see if maximum is unique")

        return(ridgeParameter)
    }
    
    ridgeParameter <- LLSolve( RTrain=RTrain,resPred=resPred,
      sdInvTrain=sdInvTrain,refs=refs, nFolds=nFolds,
		  lambdaInit=lambdaInit,iter=iter,ridgeParameter=ridgeParameter,
      LL=LL,tol=tol, max.iter=max.iter )

}

  if(ridgeParameter < tol * 10 ){
    ridgeParameter <- 0
    warning("Ridge parameter was estimated to be 0.")
  }

return(ridgeParameter)
}


################################################################################
# LLCalc: For calculating the Cross Validation penalised (log)likelihood.      #
################################################################################
LLCalc <- function(lambdas,isSing,nFolds,resPred,refs,RTrain,nVars,sdInvTrain,
    tol=1.0e-050) {

if (isSing)
  { #remove k=1 for singular cases
  lambdas<-lambdas[lambdas!=1] # will give an error, if all lambdas are =1
  }

LLPred  <- rep.int(0, times= length(lambdas))
dig <- matrix(0, nVars, nVars)
for (iLambda in 1:length(lambdas) ) {

  # Obtain the augmentation parameter for this round
	lami			<- lambdas[iLambda]
  # Find LL using K-fold validation
  for (iFold in 1:nFolds)   {
        
    isTest    <- (refs == iFold)
		l.isTest  <-	sum(isTest)
    resi      <- resPred[isTest, ,drop = FALSE]
		dig[1+0:(nVars-1)*(nVars+1)]  <- ((1-lami)/lami)				
		RHat      <- RTrain[[iFold]] +  dig
    RHatInv   <- solve( RHat ,tol=tol)
    SPredInv  <- sdInvTrain[[iFold]] %*% RHatInv %*% sdInvTrain[[iFold]]
		LLPredI   <- sum( c( resi %*% SPredInv %*% t(resi) )[1+0:(l.isTest-1)*
                  (l.isTest+1)] )
    eigVSPredInv	<- eigen(SPredInv,  only.values =TRUE)$values

    if (all(eigVSPredInv!=0) )  {
        LLPredI   <- LLPredI - sum(isTest) * sum( log( eigVSPredInv) )
    } else {
    # if there are eigenvalues=0
        eigVSPredInv[eigVSPredInv==0]<- tol
        LLPredI   <- LLPredI - l.isTest * sum( log( eigVSPredInv) )
        warning("due to matrix singularities, the penalised loglikelihood
        could not be calculated precisely, resulting in a non-precise solution" )
        }
        # LLPredI   <- LLPredI + sum(isTest) * (1-lami) * "trace"( RHatInv )
        # only include this line for penalised likelihood.

        LLPred[iLambda] <- LLPred[iLambda] + LLPredI
   }
   }
return(LLPred)

}


