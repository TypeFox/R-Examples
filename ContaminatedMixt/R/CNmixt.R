#####################################
## Bivariate Gaussian Distribution ##
#####################################

.bivnorm <- function(x,y,mu,Sigma){
  
  rho    <- Sigma[1,2]
  muX    <- mu[1]
  muY    <- mu[2]
  sigmaX <- sqrt(Sigma[1,1])
  sigmaY <- sqrt(Sigma[2,2])
  mah    <- 1/(1-rho^2)*( (x-muX)^2/sigmaX^2 - 2*rho*(x-muX)*(y-muY)/(sigmaX*sigmaY) + (y-muY)^2/sigmaY^2 )
  
  1/(2*pi*sigmaX*sigmaY*sqrt(1-rho^2))*exp(-1/2*mah)  
  
}

#################################################################
## Modified Weighted Covariance Matrix for contaminated models ##
#################################################################

.cov.wt.CN <- function(x,wt,fact){
  
  # fact contains the corrections factors due to the contamination 
  
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  else if (!is.matrix(x)) 
    stop("'x' must be a matrix or a data frame")
  if (!all(is.finite(x))) 
    stop("'x' must contain finite values only")
  p   <- ncol(x)
  n   <- nrow(x)
  mu  <- array(colSums(wt*fact/sum(wt*fact)*x),c(p),dimnames=list(paste("X.",1:p,sep="")))
  cov <- array(crossprod(sqrt(wt*fact/sum(wt))*(x-matrix(rep(mu,n),n,p,byrow=TRUE))),c(p,p),dimnames=list(paste("X.",1:p,sep=""),paste("X.",1:p,sep="")))
  return(
    list(
      center = mu,
      cov = cov      
    )
  )
  
}

################################
## M-step of the EM algorithm ##
################################

.m.step <- function(data=NULL, covtype=NULL, w=NULL, fact=matrix(1,nrow(w),ncol(w)), v=1, D=NULL, mtol=NULL, mmax=NULL) {
  
  G = ncol(w);
  d = ncol(data);
  Sk = array(0, c(d, d, G))
  gpar= list()
  for (k in 1:G) {
    gpar[[k]] = list()    
    temp = .cov.wt.CN(x=data, wt=w[,k], fact=fact[,k])
    gpar[[k]]$mu = temp$center
    if (!any(is.na( temp$cov))) 
      gpar[[k]]$sigma = temp$cov
    Sk[,,k] = temp$cov
  }
  gpar$pi = apply(w,2,mean)
  
  temp = model.type(modelname = covtype, Sk = Sk, ng = gpar$pi, D = D, mtol = mtol, mmax = mmax)
  gpar$D = temp$D
  for (k in 1:G) {
    gpar[[k]]$sigma    = temp$sigma[,,k] 
    gpar[[k]]$invSigma = temp$invSigma[,,k] 
    gpar[[k]]$logdet   = temp$logdet[k]
  }
  return(gpar)
  
}

######################################
# used in the m-step to update alpha #
######################################

.gmax <- function(alpha,j,z,v)      
  sum(z[,j]*(v[,j]*log(alpha)+(1-v[,j])*log(1-alpha)))

####################################
# used in the m-step to update eta #
####################################

.fmax <- function(eta,X,j,p,zvbad,mu,invSigma) 
  sum(-(p/2)*zvbad[,j]*log(eta)-(1/2)*zvbad[,j]*(1/eta)*mahalanobis(x=X, center=mu[,j], cov=invSigma[,,j], inverted=TRUE))

####################################################
## Mixture of Contaminated Gaussian distributions ##
## with fixed number of components                ##
####################################################

.CNmixtG <- function(
  X,			                      # matrix of data
  G,                            # number of groups
  initialization="mixt",        # initialization procedure: "kmeans", random.soft", "random.hard", "manual", or "mixt" (which is nested)
  modelname="VVV",              # one of the 14 models of Celeaux & Govaert (1995)       
  alphafix=NULL,                   # vector of dimension G with proportion of good observations in each group
  alphamin=0.5,          # vector of minimum proportions of good data 
  etafix=NULL,                     # vector of dimension G with degree of contamination in each group
  etamax=1000,                  # maximum value of eta
  seed = NULL,
  start.z = NULL,               # (n x G)-matrix of soft or hard classification: it is used only if initialization="manual"		
  start.v = NULL,               # (n x 2 x G)-array of soft or hard classification in each group: it is used only if initialization="manual"  	
  #veo = FALSE,                 # If TRUE then if the number variables in the model exceeds the number of observations the model is still fitted.
  start=0,                      # initialization for the package mixture
  ind.label=NULL,               # indexes of the labelled observations
  label=NULL,                   # groups of the labelled observations
  iter.max=1000,                # maximum number of iterations in the EM-algorithm
  threshold=1.0e-04,            # stopping rule in the Aitken rule
  eps
)
{
  
  if (is.data.frame(X)) 
    X <- as.matrix(X)
  
  n <- nrow(X)    # sample size
  p <- ncol(X)    # number of variables
  
  if (is.null(X))     stop('Hey, we need some data, please! X is null')
  if (!is.matrix(X))  stop('X needs to be in matrix form')
  if (!is.numeric(X)) stop('X is required to be numeric')
  if (n == 1)   stop('nrow(X) is equal to 1')
  if (p == 1)   stop('ncol(X) is equal to 1; This function currently only works with multivariate data (p > 1)')
  if (any(is.na(X)))  stop('No NAs allowed.')
  if (is.null(G)) stop('G is NULL')
  G <- as.integer(ceiling(G))
  if (!is.integer(G)) stop('G is not a integer')
  if (any(G < 1)) stop('G is not a positive integer')
  
  # for model-based classification
  
  lab <- NULL
  if(is.vector(label)){
    nlab           <- length(label)
    nunlab         <- n-nlab
    lab            <- numeric(n)
    lab[ind.label] <- label
  }  

  # -------------------- #
  # Number of Parameters #
  # -------------------- #
  
  if(is.null(alphafix) & is.null(etafix))
    npar <- (G-1) + p*G + .ncovpar(modelname=modelname, p=p, G=G) + 2*G
  if( (is.null(alphafix) & !is.null(etafix)) | (!is.null(alphafix) & is.null(etafix)) )
    npar <- (G-1) + p*G + .ncovpar(modelname=modelname, p=p, G=G) + G
  if(!is.null(alphafix) & !is.null(etafix))
    npar <- (G-1) + p*G + .ncovpar(modelname=modelname, p=p, G=G) 
  
  prior      <- numeric(G) # proportion of each group
  v          <- array(1,c(n,G),dimnames=list(1:n,paste("group ",1:G,sep="")))
  mu         <- array(0,c(p,G),dimnames=list(paste("X.",1:p,sep=""),paste("group ",1:G,sep="")))
  Sigma      <- array(0,c(p,p,G),dimnames=list(paste("X.",1:p,sep=""),paste("X.",1:p,sep=""),paste("group ",1:G,sep=""))) 
  invSigma   <- array(0,c(p,p,G),dimnames=list(paste("X.",1:p,sep=""),paste("X.",1:p,sep=""),paste("group ",1:G,sep=""))) 
  lambda     <- array(0,c(G),dimnames = list(paste("group ",1:G,sep="")))
  Delta      <- array(0,c(p,p,G),dimnames=list(paste("X.",1:p,sep=""),paste("X.",1:p,sep=""),paste("group ",1:G,sep="")))
  Gamma      <- array(0,c(p,p,G),dimnames=list(paste("X.",1:p,sep=""),paste("X.",1:p,sep=""),paste("group ",1:G,sep="")))
  if(is.null(alphafix))
      alpha <- rep(0.999,G) else{
      if(length(alphafix)==1) alpha <- rep(alphafix,G)
      if(length(alphafix)==G) alpha <- alphafix
      if(length(alphafix)!=1 & length(alphafix)!=G){
        alpha  <- rep(alphafix[1],G)
        warning("length of alphafix is different from G and alphafix[1] is replicated G times")
      }
    }
  if(is.null(etafix))
    eta <- rep(1.001,G) else{
      if(length(etafix)==1) eta <- rep(etafix,G)
      if(length(etafix)==G) eta <- etafix
      if(length(etafix)!=1 & length(etafix)!=G){
        eta  <- rep(etafix[1],G)
        warning("length of etafix is different from G and etafix[1] is replicated G times")
      }
    }

  if(!is.null(alphamin)){
    if(length(alphamin) == 1) alphamin <- rep(alphamin,G)
     if(length(alphamin)!=G) {
        alphamin <- rep(alphamin[1],G)
        warning("length of alphamin is different from G and alphamin[1] is replicated G times")
    }
  }
  if(!is.null(etamax)){
    if(length(etamax) == 1) etamax <- rep(etamax,G)
    if(length(etamax)!=G) {
      etamax <- rep(etamax[1],G)
      warning("length of etamax is different from G and etamax[1] is replicated G times")
    }
  }  
  
  #alpha      <- ifelse(is.null(alphafix),numeric(G),alphafix) # proportion of good observations in each group
  #eta        <- ifelse(is.null(etafix),rep(1.001,G),etafix)   # degrees of contamination
  correction <- array(0,c(n,G),dimnames=list(1:n,paste("group ",1:G,sep=""))) # factor which differentiates this model by mclust
  PXgood     <- array(0,c(n,G),dimnames=list(1:n,paste("group ",1:G,sep="")))
  PXbad      <- array(0,c(n,G),dimnames=list(1:n,paste("group ",1:G,sep="")))
  
  # ------------------------- #
  # posteriors initialization #
  # ------------------------- #
  
  if(initialization=="random.soft"){
    
    if(!is.null(seed)) 
      set.seed(seed)
    z  <- array(runif(n*G),c(n,G)) # soft posterior probabilities (no-normalized) (n x G) 
    z  <- z/rowSums(z)             # soft posterior probabilities (n x G)
    if(!is.null(start.v))    
      v <- start.v  
    
  } 
  
  if(initialization=="random.hard"){
    
    if(!is.null(seed)) 
      set.seed(seed)
    z <- t(rmultinom(n, size = 1, prob=rep(1/G,G)))  # hard posterior probabilities (n x G)
    if(!is.null(start.v))    
      v <- start.v
    
  } 
  
  if(initialization=="manual"){ 
    
    z  <- start.z      
    if(!is.null(start.v))    
      v <- start.v    
      
  }
  
  if(initialization=="kmeans"){
    
    clusters  <- kmeans(x=X, centers=G,nstart = 1)      
    z         <- mclust::unmap(clusters$cluster,1:G)
    if(!is.null(start.v))    
      v <- start.v
  } 
  
  if(initialization=="mixt"){
    
    mixturefit <- gpcm(data=X, G=G, mnames=modelname, start=start, label=lab, veo=TRUE, atol=threshold, pprogress=FALSE)
    z <- mixturefit$z
    if(!is.null(start.v))    
      v <- start.v
        
  }
  
  if(is.null(z)){
    error <- paste0("\nModel ",modelname," with G = ",G," was not estimated due to bad initialization.")
    return(list(error=error))
  }
  
  if (min(colSums(z))<=p) z <- (z+0.0000001)/rowSums(z+0.0000001)
  
  # z for labeled observations
  
  if(is.vector(label))
    z[ind.label,] <- mclust::unmap(label, G=G)
  
  # ------------ #
  # EM algorithm #
  # ------------ #
  
  # Preliminary definition of convergence criterions
  
  check     <- 0
  iteration <- 1
  loglik    <- NULL
  aloglik   <- c(0,0)
  a         <- NULL
  a         <- c(0,0)
  
  while(check<1){
    
    # ++++++ #
    # M-step #
    # ++++++ #
    
    # ------- #
    # Weights #
    # ------- #
    
    prior     <- colMeans(z)
    zv        <- z*v
    zvcompl   <- z*(1-v)
    if(is.null(alphafix)){
      if(!is.null(alphamin)){
        for(j in 1:G){
          alpha[j] <- optimize(.gmax,c(alphamin[j],1), maximum = TRUE, j=j, z=z, v=v)$maximum
        }
      }
      if(is.null(alphamin))
        alpha <- colSums(zv)/colSums(z)
    }
    priorbad <- 1 - alpha
    
    # ---------- #
    # mu & Sigma #
    # ---------- #
    
    correction <- v+(1-v)*matrix(rep(1/eta),n,G,byrow=TRUE)
    
    fitM <- .m.step(data=X, covtype=modelname, w=z, fact=correction, v=1, mtol=1e-10, mmax=10)
    for(j in 1:G){ 
      mu[,j]        <- fitM[[j]]$mu
      Sigma[,,j]    <- .fixSigma(fitM[[j]]$sigma,eps)
      invSigma[,,j] <- fitM[[j]]$invSigma
    }
    
    # ------------------- #
    # Inflation parameter #
    # ------------------- #  
    
    if(is.null(etafix)){
      for(j in 1:G)
        eta[j] <- optimize(.fmax, c(1,etamax[j]), maximum = TRUE, X=X, j=j, p=p, zvbad=zvcompl, mu=mu, invSigma=invSigma)$maximum
    }
    
    # ------- #
    # density #
    # ------- #
    
    zerostar <- .Machine$double.xmin # to avoid zero probabilities
    for(j in 1:G){
      PXgood[,j] <- (2*pi)^(-p/2)*(det(Sigma[,,j]))^(-1/2)*exp(-1/2*mahalanobis(x=X, center=mu[,j], cov=invSigma[,,j], inverted=TRUE)) 
      PXbad[,j]  <- (2*pi)^(-p/2)*(det(eta[j]*Sigma[,,j]))^(-1/2)*exp(-1/2*mahalanobis(x=X, center=mu[,j], cov=1/eta[j]*invSigma[,,j], inverted=TRUE))
      PXgood[,j] <- (PXgood[,j]<zerostar)*zerostar+(PXgood[,j]>=zerostar)*PXgood[,j]
      PXbad[,j]  <- (PXbad[,j]<zerostar)*zerostar+(PXbad[,j]>=zerostar)*PXbad[,j]
    }
    
    # ------------------------------------- # 
    # Global - Observed-data log-likelihood # 
    # ------------------------------------- #
    
    # model-based clustering
    
    if(!is.vector(label))
      llvalues <- sum(log( rowSums(matrix(rep(prior,n),n,G,byrow=TRUE)*(matrix(rep(alpha,n),n,G,byrow=TRUE)*PXgood + matrix(rep(priorbad,n),n,G,byrow=TRUE)*PXbad))))
    
    # model-based classification
    
    if(is.vector(label)){    
      llvalueslab   <- z[ind.label,]*(log(matrix(rep(prior,nlab),nlab,G,byrow=T))+log( matrix(rep(alpha,nlab),nlab,G,byrow=TRUE)*PXgood[ind.label,]+matrix(rep(priorbad,nlab),nlab,G,byrow=TRUE)*PXbad[ind.label,] ))
      llvaluesunlab <- log(rowSums(matrix(rep(prior,nunlab),nunlab,G,byrow=TRUE)*(matrix(rep(alpha,nunlab),nunlab,G,byrow=TRUE)*PXgood[-ind.label,]+matrix(rep(priorbad,nunlab),nunlab,G,byrow=TRUE)*PXbad[-ind.label,])))
      llvalues      <- sum(llvalueslab)+sum(llvaluesunlab)    
    }
    
    loglik[iteration] <- llvalues
    
    # --------------------------- #
    # Aitken's Stopping Criterion #
    # --------------------------- #
    
    if(iteration>2 & G > 1){
      if(abs(loglik[iteration-1]-loglik[iteration-2])>0){
        a[iteration-1]      <- (loglik[iteration]-loglik[iteration-1])/(loglik[iteration-1]-loglik[iteration-2])
        aloglik[iteration]  <- loglik[iteration-1]+(1/(1-a[iteration-1])*(loglik[iteration]-loglik[iteration-1]))
        if(abs(aloglik[iteration]-loglik[iteration])<threshold) 
          check <- 1
      }
      else
        check <- 1    
    }
    
    if(iteration==iter.max | G==1) 
      check <- 1
    
    cat("*")
    iteration <- iteration + 1
    
    # ++++++ #
    # E-Step #
    # ++++++ #
    
    z.num  <- matrix(rep(prior,n),n,G,byrow=TRUE)*(matrix(rep(alpha,n),n,G,byrow=TRUE)*PXgood+matrix(rep(priorbad,n),n,G,byrow=TRUE)*PXbad)  # (n x G)
    z.den  <- rowSums(z.num)                        # n-vector
    z      <- z.num/matrix(rep(z.den,G),ncol=G)     # (n x G)
    
    # z adjustment
    if (min(colSums(z))<=p) z <- (z+0.0000001)/rowSums(z+0.0000001)
    
    if(is.vector(label))
      z[ind.label,] <- mclust::unmap(label,G=G)
    
    numgood <- matrix(rep(alpha,n),n,G,byrow=TRUE)*PXgood
    v.den   <- matrix(rep(alpha,n),n,G,byrow=TRUE)*PXgood + matrix(rep(priorbad,n),n,G,byrow=TRUE)*PXbad
    v       <- numgood/v.den
    
  }
  
  cat("\n")
  finalloglik <- loglik[iteration-1] 
  
  # **************************** #
  # The EM-algorithm is finished #
  # **************************** #
  
  # --------------------- #
  # Classification Matrix #
  # --------------------- #
  
  group <- apply(z,1,which.max)
  innergroup  <- numeric(n)
  for(i in 1:n)
    innergroup[i] <- ifelse(v[i,group[i]]<0.5,"bad","*")
  detection <- data.frame(group=group,innergroup=innergroup)
  
  # -------------------- #
  # Information Criteria #
  # -------------------- #
  
  df <- npar
  IC <- list()
  IC$AIC   <- 2*finalloglik - df*2
  IC$BIC   <- 2*finalloglik - df*log(n)
  IC$AIC3  <- 2*finalloglik - df*3  
  IC$AICc  <- IC$AIC - (2*df*(df+1))/(n-df-1)
  IC$AICu  <- ifelse(n/(n-df-1)>0,IC$AICc - n*log(n/(n-df-1)),NA)
  IC$CAIC  <- 2*finalloglik - df*(1+log(n))
  IC$AWE   <- 2*finalloglik - 2*df*(3/2+log(n))  
  z.const  <- (z<.Machine$double.xmin)*.Machine$double.xmin+(z>.Machine$double.xmin)*z   # vincolo per evitare i NaN nel calcolo di tau*log(tau)
  hard.z   <- (matrix(rep(apply(z,1,max),G),n,G,byrow=F)==z)*1
  ECM      <- ifelse(is.vector(label),sum(hard.z[-ind.label,]*log(z.const[-ind.label,])),sum(hard.z*log(z.const)))
  IC$ICL   <- IC$BIC+ECM
  
  # ------- #
  # results #
  # ------- #
  
  result <- list(
    model  = modelname,
    npar      = npar,
    X         = X,            
    G         = G,            
    p         = p,            
    n         = n,            
    prior     = prior,
    alpha     = alpha,
    mu        = mu,
    Sigma     = Sigma,
    eta       = eta,
    iter.stop = iteration,
    posterior = z,
    v         = v,
    ind.label = ind.label,               
    label     = label,                   
    group     = group,
    detection = detection,
    loglik    = finalloglik,
    IC        = IC,
    call      = match.call()
  )
  
  class(result) <- "CNmixt"
  return(result)
}

.fixSigma <- function(sigma,eps){
  es <- eigen(sigma)
  es$values[es$values<eps] <- eps
  es$vectors %*% diag(es$values) %*% t(es$vectors)
}
