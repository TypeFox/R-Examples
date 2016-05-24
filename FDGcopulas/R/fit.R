## Empirical pairwise dependence coefficients 
# (matrix)
corFDG <- function(x, depcoefType){
  stopifnot(depcoefType %in% c("pearson","spearman","kendall","utdc"))
  d <- ncol(x); p <- d*(d-1)/2
  if(depcoefType=="utdc"){
    utdc <- matrix(ncol=d,nrow=d)
    for(i in 1:(d-1)){
      for(j in (i+1):d){
        utdc[j,i]<-utdc[i,j]<-3-1/( 1 - mean( apply(x[,c(i,j)], 1, max) ) )
      }
    }
    diag(utdc) <- rep(1,d)
    utdc
  }else{
      cor(x, method=depcoefType)
  }
}






## Jacobian matrix of the map that to the parameter vector associates
# the vector whose coordinates are the pairwise dependence
# coefficients
# - theta : multivariate parameter vector
# - depcoefType : type of dependence coefficient defining the above
# mentionned map
JacobianFDG <- function(copula, depcoefType){
  theta <- copula@parameters
  tfun <- dispatch(copula, depcoefType)$tfun
  d <- length(theta); p <- d*(d-1)/2
  Jacob <- matrix(nrow=p, ncol=d, 0)
  k <- 0
  for(i in 1:(d-1)){
    for(j in (i+1):d){
      k <- k+1
      partial.theta.i <- tfun(theta[i],theta[j]) # deriv w.r.t. theta_i
      partial.theta.j <- tfun(theta[j],theta[i]) # deriv w.r.t. theta_j
      Jacob[k,c(i,j)] <- c(partial.theta.i,partial.theta.j)
    }
  }
  Jacob
}



## Loss function,
# to be minimized when estimating the parameters. Returns a list whose
# components are loss, the actual loss, and gr, its gradient w.r.t.
# the parameters
# - 'dcData' is the matrix of pairwide empirical dependence coefficients.
# - depcoefType : type of dependence coefficient 
# - 'W' is the weight matrix.
lossFDG <- function(copula, dcData, W, depcoefType){
    theta <- copula@parameters
    d <- length(theta)
    p <- d*(d-1)/2
    ## rhoFoo <- function(a,b){
    ##   rhoFDG(FDGcopula(copula@family,c(a,b),copula@extremevalue))[1,2]
    ## }
    ## tauFoo <- function(a,b){
    ##   tauFDG(FDGcopula(copula@family,c(a,b),copula@extremevalue))[1,2]
    ## }
    ## utdcFoo <- function(a,b){
    ##   utdcFDG(FDGcopula(copula@family,c(a,b),copula@extremevalue))[1,2]
    ## }
    Foo <- switch(depcoefType,
                  "spearman"=dispatch(copula, depcoefType)$rho,
                  "kendall"=dispatch(copula, depcoefType)$tau,
                  "utdc"=dispatch(copula, depcoefType)$utdc)
    DiffMat <- dcData-outer(theta,theta,Vectorize(Foo))
    DiffVec <- numeric(length=p)
    k <- 0
    for(i in 1:(d-1)){
        for(j in (i+1):d){
            k <- k+1
            DiffVec[k] <- DiffMat[i,j]
        }
    }
    loss <- t(DiffVec)%*%W%*%DiffVec
    Jacob <- JacobianFDG(copula, depcoefType)
    gr <- -2*t(DiffVec)%*%W%*%Jacob
    list(loss=loss, gr=gr)
}



## Weighted least squares inference based on dependence coefficients.
# - 'dcData' is the matrix of empirical dependence coefficients,
# - 'loss' is the loss function,
# The global optimization problem is solved by brute force:
# we provide 'N' initialisations to 'optim' and see which gives the
# best results. The objective values of each 'optim' call are stored
# in the vector 'lossImprove'. The global minimizer is stored in 'para'.
wlsdc <- function(copula, dcData,  depcoefType, nbInit, W, method){
  d <- length(copula@parameters)
  objectiveValue = 0
  optimPara = rep(0, d)
  improvementOV = 0
  sob <- sobol(nbInit,dim=d)
  upper <- copula@parameterrange[2]
  lower <- copula@parameterrange[1]

  loss <- function(theta){
    lossFDG(FDGcopula(copula@family, theta, copula@extremevalue, 
                      checkbounds=FALSE),
            dcData, W, depcoefType)$loss
  }
  gr <- function(theta){
    lossFDG(FDGcopula(copula@family, theta, copula@extremevalue,
                      checkbounds=FALSE),
            dcData, W, depcoefType)$gr
  }

  convergenceDiagnostic <- list()
  for(k in 1:nbInit){
    optimOut <- if(method=="L-BFGS-B") {
      optim(sob[k,]*upper, loss, gr=gr, method="L-BFGS-B",
                     lower=lower, upper=upper)
    }else{
      optim(sob[k,]*upper, loss, gr=gr, method=method)
    }

    optimPara <- rbind(optimPara, optimOut$par)
    objectiveValue <- c(objectiveValue, optimOut$value)
    improvementOV <- c(improvementOV, min(objectiveValue[-1]))
    convergenceDiagnostic[[k]] <- list(counts=optimOut$counts,
                                       convergence=optimOut$convergence,
                                       message=optimOut$message)
  }

  objectiveValue <- objectiveValue[-1]
  improvementOV <- improvementOV[-1]
  optimPara <- optimPara[-1,]
  if(is.null(dim(optimPara))){
    para <- optimPara
  }else{
    pos <- which(objectiveValue == min(objectiveValue))
    if(length(pos)>1) pos <- pos[1] # (to return a single estimate)
    para <- optimPara[pos,]
  }

  list(estimate=para, optimalValues=improvementOV,
       convergenceDiagnostic=convergenceDiagnostic)
}



## Asymptotic (co)variance of dependence coefficients #
# Asymptotic variance-covariance matrix \Sigma of the vector whose coordinates are the considered dependence coefficients. The calculation is done by simulating 'nb.rep' datasets of size 'nb.obs' according to the copula 'copula'. The implemented method are different wether tail coefficients are used or not.
asymp.vcov <- function(nb.rep, nb.obs, copula, depcoefType,
                       sizeSubSample=10000){
  theta <- copula@parameters
  if(copula@extremevalue){
    estimCov_ev_CPP(rFDG(nb.rep, copula, sizeSubSample))
  }else{    
    d <- length(theta); p <- d*(d-1)/2
    depc <- matrix(nrow=nb.rep, ncol=p)
    for(k in 1:nb.rep){
      rep <- rFDG(nb.obs, copula)
      P <- corFDG(rep, depcoefType)
      depc[k,] <- P[lower.tri(P)]
    }
    nb.obs * cov(depc)
  }
}



## Fit :
# - 'data' is the data set
# - 'copula' is the object of class 'FDGcopula'
# - 'depcoefType' is the chosen type of dependence coefficient,
# currently one of 'spearman', 'kendall', or 'utdc' (the upper
# tail dependence coefficient)
# - 'W' is the weight matrix (simply 'NA' if no weights are wanted)
# - 'method' is the 'method' argument passed to 'stats::optim'
# - 'nbInit' is the number of initialisations wanted to solve the
# global optimization problem.
# - dcData is the matrix consisting of the pairwise empirical
# dependence coefficients. If NA, it is calculated according to the
# value of depcoefType.

fitFDG <- function(FDGcopula, data, depcoefType="spearman", nbInit=1,
                   W=NA, method="L-BFGS-B", estimate.variance = TRUE,
                   nb.rep=100, nb.obs=100, dcData=NA, sizeSubSample=10000){
  if(!depcoefType%in%c("spearman","kendall","utdc")){
    stop("'depcoefType' should be one of 'spearman', 'kendall', or 'utdc'")
  }else{}
  if(FDGcopula@extremevalue==TRUE & FDGcopula@family=="exponential"){
    stop("The extreme-value copula coming from the FDG copula with exponential generators is the INDEPENDENCE copula!")
  }else{}
  
  d <- length(FDGcopula@parameters)
  p <- d*(d-1)/2
  if(is.na(W)==TRUE){
    W <- diag(1,p)
  }else{}

  if(is.na(sum(dcData))){
    dcData <- corFDG(data, depcoefType)
    if(depcoefType=="utdc" & FDGcopula@extremevalue==FALSE){
      stop("The upper tail dependence coefficient estimator implemented by default is theoretically well grounded ONLY for extreme-value copulas. If you still wanna do this way, use your own estimator by providing 'dcData'")
    }else{}
  }else{}
  
  optimResults <- wlsdc(FDGcopula, dcData, depcoefType, nbInit, W, method)
  Sigma <- if(estimate.variance){
      asymp.vcov(nb.rep, nb.obs, FDGcopula, depcoefType, sizeSubSample)
  }else{
      matrix(ncol=p,nrow=p,rep(NA,p*p))
  }
  fittedFDG <- FDGcopula(FDGcopula@family, optimResults$estimate,
                         FDGcopula@extremevalue)
  J <- JacobianFDG(fittedFDG, depcoefType)
  invJtJ <- try( solve(t(J)%*%J,diag(1,d)) )
  Xi <- if(class(invJtJ)=="try-error"){
    matrix(ncol=d,nrow=d,NA)
  }else{
    invJtJ%*%t(J)%*%Sigma%*%J%*%invJtJ
  }
  new(Class="fitFDG", estimate=optimResults$estimate,
      var.est=Xi, optimalvalues=optimResults$optimalValues,
      convergence=optimResults$convergenceDiagnostic,
      FDGcopula=fittedFDG)
}
