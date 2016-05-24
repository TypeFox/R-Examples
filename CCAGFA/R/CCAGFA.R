#
# Implementation of Bayesian models for canonical correlation analysis (CCA),
# inter-battery factor analysis (BIBFA) and group factor analysis (GFA).
#

GFAexperiment <- function(Y,K,opts,Nrep=10) {
  #
  # A wrapper for running the GFA model Nrep times
  # and choosing the final model based on the best
  # lower bound. This is the recommended way of applying
  # the algorithm.
  #
  # See GFA() for description of the inupts.
  #
  if(Nrep==1) {
    return(GFA(Y,K,opts))
  }
  
  lb <- vector() # Lower bounds
  models <- vector("list")
  for(rep in 1:Nrep){
    model <- GFA(Y,K,opts)
    models[[rep]] <- model
    lb[rep] <- tail(model$cost,1)
    if(opts$verbose==1) {
      print(paste("Run ",rep,"/",Nrep,": ",length(model$cost)," iterations",
                  " with final cost ",lb[rep],sep=""))
    }
    }
  
  k <- which.max(lb)
  model <- models[[k]]
  
  return(model)
  }

CCAexperiment <- function(Y,K,opts,Nrep=10) {
  #
  # A wrapper for running the CCA/BIBFA model Nrep times
  # and choosing the final model based on the best
  # lower bound. This is the recommended way of applying
  # the algorithm.
  #
  # See CCA() for description of the inputs.
  #
  
  if(length(Y)!=2) {
    print("Canonical correlation analysis (CCA) and Bayesian inter-battery
          factor analysis (BIBFA) can only be computed for two co-occurring
          data matrices. For multiple co-occurring data sources, use GFA()
          and GFAexperiment() instead.")
  } else {
    return(GFAexperiment(Y,K,opts,Nrep))
  }
}

CCA <- function(Y,K,opts) {
  #
  # A wrapper function for computing Bayesian CCA and BIBFA.
  # Calls the same function that is used for GFA.
  #
  if(length(Y)!=2) {
    print("Canonical correlation analysis (CCA) and Bayesian inter-battery
          factor analysis (BIBFA) can only be computed for two co-occurring
          data matrices. For multiple co-occurring data sources, use GFA()
          and GFAexperiment() instead.")
  } else {
    return(GFA(Y,K,opts))
  }
}

GFA <- function(Y,K,opts) {
  #
  # The main function for Bayesian group factor analysis
  #
  # Inputs:
  #   Y    : List of M data matrices. Y[[m]] is a matrix with
  #          N rows (samples) and D_m columns (features). The
  #          samples need to be co-occurring.
  #          NOTE: All of these should be centered, so that the mean
  #                of each feature is zero
  #          NOTE: The algorithm is roughly invariant to the scale
  #                of the data, but extreme values should be avoided.
  #                Data with roughly unit variance or similar scale
  #                is recommended.
  #   K    : The number of components
  #   opts : List of options (see function getDefaultOpts())
  #
  # Output:
  # The trained model, which is a list that contains the following elements:
  #   Z    : The mean of the latent variables; N times K matrix
  #   covZ : The covariance of the latent variables; K times K matrix
  #   ZZ   : The second moments ZZ^T; K times K matrix
  #
  #   W    : List of the mean projections; D_i times K matrices
  #   covW : List of the covariances of the projections; D_i times D_i matrices
  #   WW   : List of the second moments WW^T; K times K matrices
  #
  #   tau  : The mean precisions (inverse variance, so 1/tau gives the
  #          variances denoted by sigma in the paper); M-element vector
  #
  #   alpha: The mean precisions of the projection weights, the
  #          variances of the ARD prior; M times K matrix
  #
  #   U,V,u.mu,v.mu: The low-rank factorization of alpha.
  #
  #   cost : Vector collecting the variational lower bounds for each
  #          iteration
  #   D    : Data dimensionalities; M-element vector
  #   datavar   : The total variance in the data sets, needed for
  #               GFAtrim()
  #   addednoise: The level of extra noise as in opts$addednoise
  #
  
  # Check that data is centered
  if(!all(unlist(lapply(Y,colMeans)) < 1e-7) & opts$verbose == 2) {
    print("Warning: The data does not have zero mean.")
  }
  # Check the number of views
  if(length(Y)==1 & opts$verbose == 2){
    print("Warning: The number of data sets must be larger than 1.")
  }

  #
  # Store dimensionalities of data sets 
  #
  M <- length(Y)              # The number of views; M=2 corresponds to BIBFA
  D <- unlist(lapply(Y,ncol)) # Collect the number of features in vector D
  Ds <- sum(D)                # Total number of features
  N <- nrow(Y[[1]])           # The number of samples
  datavar <- vector()                 # The total variance of the data, needed
  for(m in 1:M)                       #     for scaling in the initialization
    datavar[m]=sum(apply(Y[[m]],2,var)) #     and for GFAtrim()

  R <- opts$R                 # Model rank R: either "full" or an integer from 0 to min(M,K)
  if(R >= min(c(M,K)) && is.double(R)) {  # if R >= min(M,K) convert to full
    if(opts$verbose==2)
      print("The rank corresponds to full rank solution.")
    R <- "full"
  }
  if(R != "full") {
    if(opts$verbose==2)
      print("NOTE: Optimization of the rotation is not supported for low rank model.")
    opts$rotate <- FALSE
  }
  
  # Some constants for speeding up the computation
  const <- - N*Ds/2*log(2*pi) # Constant factors for the lower bound
  Yconst <- unlist(lapply(Y,function(x){sum(x^2)}))
  id <- rep(1,K)              # Vector of ones for fast matrix calculations
  alpha_0 <- opts$prior.alpha_0   # Easier access for hyperprior values
  beta_0 <- opts$prior.beta_0
  alpha_0t <- opts$prior.alpha_0t
  beta_0t <- opts$prior.beta_0t
  
  #
  # Initialize the model randomly; other initializations could
  # be done, but overdispersed random initialization is quite good.
  # 
  
  # Latent variables
  Z <- matrix(rnorm(N*K,0,1),N,K) # The mean 
  covZ <- diag(1,K)               # The covariance
  ZZ <- covZ + covZ*N             # The second moments
  
  # ARD and noise parameters
  alpha <- matrix(1,M,K)          # The mean of the ARD precisions
  logalpha <- matrix(1,M,K)       # The mean of <\log alpha >
  if(R=="full") {
    b_ard <- matrix(1,M,K)        # The parameters of the Gamma distribution
    a_ard <- alpha_0 + D/2        #     for ARD precisions
    digammaa_ard <- digamma(a_ard)  # Constant for the lower bound
  }
  tau <- rep(opts$init.tau,M)     # The mean noise precisions
  a_tau <- alpha_0t + N*D/2       # The parameters of the Gamma distribution
  b_tau <- rep(0,M)               #     for the noise precisions
  digammaa_tau <- digamma(a_tau)  # Constants needed for computing the lower bound
  lgammaa_tau <- -sum(lgamma(a_tau))
  lb.pt.const <- -M*lgamma(alpha_0t) + M*alpha_0t*log(beta_0t)
  
  # Alpha needs to be initialized to match the data scale
  for(m in 1:M) {
    alpha[m,] <- K*D[m]/(datavar[m]-1/tau[m])
  }
  
  # The projections
  # No need to initialize the projections randomly, since their updating
  # step is the first one; just define the variables here
  W <- vector("list",length=M)    # The means
  if(!opts$low.mem)
    covW <- vector("list",length=M) # The covariances
  else
    covW <- diag(1,K)
  
  WW <- vector("list",length=M)   # The second moments
  for(m in 1:M) {
    W[[m]] <- matrix(0,D[m],K)
    if(!opts$low.mem){
      covW[[m]] <- diag(1,K)
      WW[[m]] <- crossprod(W[[m]]) + covW[[m]]*D[m]
    } else {
      WW[[m]] <- crossprod(W[[m]]) + covW*D[m]
    }
  }
  
  # Rotation parameters (full-rank only)
  if(opts$rotate) {
    Rot <- diag(K)      # The rotation matrix
    RotInv <- diag(K)   # Its inverse
    r <- as.vector(Rot) # Vectorized version of R
    
    # parameter list for the optimization function (see ?optim)
    par <- list(K=K,D=D,Ds=Ds,N=N,WW=WW,ZZ=ZZ,M=M)
  }

  # Use R-rank factorization of alpha.
  if(R!="full") {
    U <- matrix(abs(rnorm(M*R)),M,R)
    lu <- length(U)
    u.mu <- rep(0,M)
    V <- matrix(abs((rnorm(K*R))),K,R)
    lv <- length(V)
    v.mu <- rep(0,K)
    
    lambda <- opts$lambda
    
    x <- c(as.vector(U),as.vector(V),u.mu,v.mu)
    x <- rnorm(length(x))/100 
    
    par.uv <- list(getu=1:lu,getv=(lu+1):(lu+lv),getumean=(lu+lv+1):(lu+lv+M),
                   getvmean=(lu+lv+M+1):length(x),M=M,K=K,R=R,D=D,lambda=lambda)
    par.uv$w2 <- matrix(0,M,K)
  }
  
  cost <- vector()  # For storing the lower bounds
  
  #
  # The main loop
  #
  for(iter in 1:opts$iter.max) {
    
    # Check if some components need to be removed
    keep <- which(colMeans(Z^2) > 1e-7)
    if(length(keep)!=K && opts$dropK) {
      K <- length(keep)
      if(K==0)
        stop("Shut down all components, no structure found in the data.")
      id <- rep(1,K)
      Z <- Z[,keep,drop=F]
      covZ <- covZ[keep,keep,drop=F]
      ZZ <- ZZ[keep,keep,drop=F]
      for(m in 1:M) {
        W[[m]] <- W[[m]][,keep,drop=FALSE]
        if(!opts$low.mem)
          covW[[m]] <- covW[[m]][keep,keep,drop=F]
        WW[[m]] <- WW[[m]][keep,keep,drop=F]
      }
      
      alpha <- alpha[,keep,drop=F]
      logalpha <- logalpha[,keep,drop=F]
      
      if(R!="full") {
        V <- V[keep,,drop=F]
        v.mu <- v.mu[keep]
        x <- c(as.vector(U),as.vector(V),u.mu,v.mu)
        lv <- length(V)
        par.uv$K <- K
        par.uv$getv<- (lu+1):(lu+lv)
        par.uv$getumean <- (lu+lv+1):(lu+lv+M)
        par.uv$getvmean <- (lu+lv+M+1):length(x)
        par.uv$w2 <- matrix(0,M,K)
      } else {
        b_ard <- matrix(1,M,K) # The parameters of the Gamma distribution
      }
      if(opts$rotate) {
        par$K <- K
      }
    }
    #
    # Update the projections
    #
    lb.qw <- rep(NA,M) # Computes also the determinant of covW needed for the
                       # lower bound
    for(m in 1:M) {
      # Efficient and robust way of computing
      # solve(diag(alpha) + tau * ZZ^T)
      tmp <- 1/sqrt(alpha[m,])

      cho <- chol(outer(tmp,tmp)*ZZ + diag(1/tau[m],K))
      det <- -2*sum(log(diag(cho))) - sum(log(alpha[m,])) - K*log(tau[m])
      lb.qw[m] <- det
      if(!opts$low.mem) {
        covW[[m]] <- 1/tau[m] * outer(tmp,tmp) * chol2inv(cho)
        W[[m]] <- crossprod(Y[[m]],Z)%*%covW[[m]]*tau[m]
        WW[[m]] <- crossprod(W[[m]]) + covW[[m]]*D[m]
      } else {
        covW <- 1/tau[m] * outer(tmp,tmp) * chol2inv(cho)
        W[[m]] <- crossprod(Y[[m]],Z)%*%covW*tau[m]
        WW[[m]] <- crossprod(W[[m]]) + covW*D[m]
      }
    }
    
    # 
    # Update the latent variables
    #
    
    # Efficient and robust way of computing
    # solve(diag(1,K) + tau * WW^t)
    covZ <- diag(1,K)
    for(m in 1:M) {
      covZ <- covZ + tau[m]*WW[[m]]
    }
    cho <- chol(covZ)
    covZ <- chol2inv(cho)
    det <- -2*sum(log(diag(cho)))
    lb.qx <- det
    
    Z <- Z*0
    for(m in 1:M)
      Z <- Z + Y[[m]]%*%W[[m]]*tau[m]
    Z <- Z%*%covZ
    ZZ <- crossprod(Z) + N*covZ
    
    #
    # Optimization of the rotation (only start after the first
    # iteration)
    #
    if(R=="full" & opts$rotate & iter > 1) {
      # Update the parameter list for the optimizer
      par$WW <- WW
      par$ZZ <- ZZ
      
      # Always start from the identity matrix, i.e. no rotation
      r <- as.vector(diag(K))
      if(opts$opt.method == "BFGS") {
        r.opt <- try(optim(r,E,gradE,par,method="BFGS",
                           control=list(reltol=opts$opt.bfgs.crit,
                                        maxit=opts$opt.iter)), silent=TRUE)
      }
      if(opts$opt.method== "L-BFGS") {
        r.opt <- try(optim(r,E,gradE,par,method="L-BFGS-B",
                           control=list(maxit=opts$opt.iter,
                                        factr=opts$lbfgs.factr)), silent=TRUE)
      }
      
      # For large K the optimizer occasionally fails due to problems
      # in the svd() routine
      if(inherits(r.opt,"try-error")) {
        print("Failure in optimizing the rotation. Turning the rotation off.")
        opts$rotate <- FALSE
      }else{
        # Update the parameters involved in the rotation: 
        Rot <- matrix(r.opt$par,K)
        eS <- svd(Rot)
        det <- sum(log(eS$d))
        RotInv <- tcrossprod( eS$v*outer(id, 1/eS$d), eS$u)
        
        Z <- tcrossprod(Z,RotInv)
        covZ <- tcrossprod( RotInv%*%covZ, RotInv)
        ZZ <- crossprod(Z) + N*covZ
        
        lb.qx <- lb.qx - 2*det
        
        for(m in 1:M) {
          if(!opts$low.mem) {
            W[[m]] <- W[[m]]%*%Rot
            covW[[m]] <- crossprod(Rot,covW[[m]])%*%Rot
            WW[[m]] <- crossprod(W[[m]]) + covW[[m]]*D[m]
          } else {
            #covW[[m]] is not stored, so it needs to be computed before rotation
            covW <- (WW[[m]] - crossprod(W[[m]]))/D[m]
            W[[m]] <- W[[m]]%*%Rot
            covW <- crossprod(Rot,covW)%*%Rot
            WW[[m]] <- crossprod(W[[m]]) + covW*D[m]
          }
          
          lb.qw[m] <- lb.qw[m] + 2*det
        }
      }
    }
    
    #
    # Update alpha, the ARD parameters
    #
    if(R == "full"){
      for(m in 1:M){
        tmp <- beta_0 + diag(WW[[m]])/2
        alpha[m,] <- a_ard[m]/tmp
        b_ard[m,] <- tmp
      }
    }else{
      for(m in 1:M)
        par.uv$w2[m,] <- diag(WW[[m]])

      res <- try(optim(x,fn=Euv,gr=gradEuv,par.uv,method="L-BFGS-B",
                            lower=c(rep(-sqrt(500/R),M*R+K*R),rep(-50,M+K)),
                            upper=c(rep(sqrt(500/R),M*R+K*R), rep(50,M+K)),  
                            control=list(fnscale=-1,maxit=opts$opt.iter,
                                         factr=opts$lbfgs.factr)),silent=TRUE)
      if(inherits(res,"try-error")){
        cost[iter] <- NA
        stop("Problems in optimization. Try a new initialization.")
        # terminate the algorithm (next model to learn) 
      }
      x <- res$par
      U <- matrix(res$par[par.uv$getu],par.uv$M,par.uv$R)
      V <- matrix(res$par[par.uv$getv],par.uv$K,par.uv$R)
      u.mu <- res$par[par.uv$getumean]
      v.mu <- res$par[par.uv$getvmean]
      alpha <- exp(tcrossprod(U,V) + outer(u.mu,rep(1,K)) +outer(rep(1,M),v.mu))
      
    }
    
    #
    # Update tau, the noise precisions
    #
    for(m in 1:M) {
      b_tau[m] <- beta_0t + ( Yconst[m] + sum(WW[[m]]*ZZ)
                              - 2*sum(Z*(Y[[m]]%*%W[[m]])))/2
    }
    tau <- a_tau/b_tau
    
    #
    # Calculate the lower bound.
    # Consists of calculating the likelihood term and 
    # KL-divergences between the factorization and the priors
    #
    
    # The likelihood
    logtau <- digammaa_tau - log(b_tau)
    if(R == "full") {
      for(m in 1:M) {
        logalpha[m,] <- digammaa_ard[m] - log(b_ard[m,])
      }
    } else {
      logalpha <- log(alpha)
    }
    
    lb.p <- const + N*crossprod(D,logtau)/2 - crossprod(b_tau - beta_0t,tau)
    lb <- lb.p
    
    # E[ ln p(Z)] - E[ ln q(Z) ]
    lb.px <- - sum(diag(ZZ))/2
    lb.qx <- - N*lb.qx/2 - N*K/2 
    lb <- lb + lb.px - lb.qx
    
    # E[ ln p(W)] - E[ ln q(W)]
    if(R == "full") {
      lb.pw <- 0
      for(m in 1:M) {
        lb.pw <- lb.pw +
          D[m]/2*sum(logalpha[m,]) -  sum(diag(WW[[m]])*alpha[m,])/2
      }
    } else {
      lb.pw <- res$value
    }
    for(m in 1:M)
      lb.qw[m] <- - D[m]*lb.qw[m]/2 -D[m]*K/2
    lb <- lb + lb.pw - sum(lb.qw)
    
    # E[ln p(alpha)] - E[ln q(alpha)]
    if(R == "full") {
      lb.pa <- M*K*( -lgamma(alpha_0) + alpha_0*log(beta_0) ) +
        (alpha_0-1)*sum(logalpha) - beta_0*sum(alpha)
      lb.qa <- -K*sum(lgamma(a_ard)) + sum(a_ard*rowSums( log(b_ard) )) +
        sum((a_ard-1)*rowSums(logalpha)) - sum(b_ard*alpha)
      lb <- lb + lb.pa - lb.qa
    }
    
    # E[ln p(tau)] - E[ln q(tau)]
    lb.pt <- lb.pt.const + sum((alpha_0t-1)*logtau) - sum(beta_0t*tau)
    lb.qt <- lgammaa_tau + crossprod(a_tau,log(b_tau)) +
      crossprod((a_tau-1),logtau) - crossprod(b_tau,tau)
    lb <- lb + lb.pt - lb.qt
    
    # Store the cost function
    cost[iter] <- lb
    
    if(opts$verbose==2) {
      print(paste("Iteration:",iter,"/ cost:",lb,"/ K",K))
    }
    # Convergence if the relative change in cost is small enough
    if(iter>1){
      diff <- cost[iter]-cost[iter-1]
      if(abs(diff)/abs(cost[iter]) < opts$iter.crit | iter == opts$iter.max) {
        break
      }
    }
  } # the main loop of the algorithm ends
  
  # Add a tiny amount of noise on top of the latent variables,
  # to supress possible artificial structure in components that
  # have effectively been turned off
  Z <- Z + opts$addednoise*matrix(rnorm(N*K,0,1),N,K) %*% chol(covZ)
  
  # return the output of the model as a list
  if(R=="full") {
    list(W=W,covW=covW,ZZ=ZZ,WW=WW,Z=Z,covZ=covZ,tau=tau,alpha=alpha,cost=cost,
         D=D,K=K,addednoise=opts$addednoise,datavar=datavar,iter=iter,R=R)
  } else {
    list(W=W,covW=covW,ZZ=ZZ,WW=WW,Z=Z,covZ=covZ,tau=tau,alpha=alpha,cost=cost,
         D=D,K=K,addednoise=opts$addednoise,datavar=datavar,iter=iter,R=R,
         U=U,V=V,u.mu=u.mu,v.mu=v.mu)
  }
}

E <- function(r,par) {
  #
  # Evaluates the (negative) cost function value wrt the transformation
  # matrix R used in the generic optimization routine
  #
  R <- matrix(r,par$K)
  eS <- svd(R)
  
  val <- -sum(par$ZZ*tcrossprod(eS$u*outer(rep(1,par$K),1/eS$d)))/2
  val <- val + (par$Ds-par$N)*sum(log(eS$d))
  for(m in 1:par$M) {
    val <- val - par$D[m]*sum( log( colSums(R*(par$WW[[m]]%*%R)) ))/2
  }
  return(-val)
}

gradE <- function(r,par) {
  #
  # Evaluates the (negative) gradient of the cost function E()
  #
  R <- matrix(r,par$K)
  eS <- svd(R)
  Rinv <- tcrossprod(eS$v*outer(rep(1,par$K),1/eS$d), eS$u)
  gr <- as.vector(tcrossprod(tcrossprod(eS$u*outer(rep(1,par$K),1/eS$d^2),eS$u)%*%
                               par$ZZ + diag(par$Ds - par$N,par$K), Rinv) )
  
  tmp1 <- par$WW[[1]]%*%R
  tmp2 <- 1/colSums(R*tmp1)
  tmp1 <- par$D[1]*as.vector( tmp1*outer(rep(1,par$K),tmp2) )
  gr <- gr - tmp1
  for(m in 2:par$M){
    tmp1 <- par$WW[[m]]%*%R
    tmp2 <- 1/colSums(R*tmp1)
    tmp1 <- par$D[m]*as.vector( tmp1*outer(rep(1,par$K),tmp2) )
    gr <- gr - tmp1
  }
  return(-gr)
}

Euv <- function(x,par){
  #
  # Evaluates the cost function value wrt the low-rank
  # factorization of alpha used in the generic optimization routine
  #
  U <- matrix(x[par$getu],par$M,par$R)
  V <- matrix(x[par$getv],par$K,par$R)
  u.mu <- x[par$getumean]
  v.mu <- x[par$getvmean]
  logalpha <- tcrossprod(U,V) + outer(u.mu,rep(1,par$K)) + outer(rep(1,par$M),v.mu)
  
  E <- sum(crossprod(par$D,logalpha))-sum(par$w2*exp(logalpha))
  if(par$lambda!=0)
    E <- E - par$lambda*(sum(V^2) + sum(U^2))
  return(E/2)
}

gradEuv <- function(x,par){
  #
  # Evaluates the gradient of the cost function Euv()
  #
  U <- matrix(x[par$getu],par$M,par$R)
  V <- matrix(x[par$getv],par$K,par$R)
  u.mu <- x[par$getumean]
  v.mu <- x[par$getvmean]
  alphaiAlphaw2 <- outer(par$D,rep(1,par$K)) - 
    exp(tcrossprod(U,V) + outer(u.mu,rep(1,par$K)) + outer(rep(1,par$M),v.mu))*par$w2
  
  gradU <- alphaiAlphaw2%*%V
  gradV <- crossprod(alphaiAlphaw2,U)
  if(par$lambda!=0) {
    gradU <- gradU - par$lambda*2*U
    gradV <- gradV - par$lambda*2*V
  }
  grad.umean <- rowSums(alphaiAlphaw2)
  grad.vmean <- colSums(alphaiAlphaw2)
  grad <- c(as.vector(gradU),as.vector(gradV),grad.umean,grad.vmean)
  return(grad/2)
}

getDefaultOpts <- function(){
  #
  # A function for generating a default set of parameters.
  #
  # To run the algorithm with other values:
  #   opts <- getDefaultOpts()
  #   opts$opt.method <- "BFGS"
  #   model <- GFA(Y,K,opts)
  
  #
  # Rank of GFA's hierarhical low-rank ARD prior. Possible values are "full" 
  # and all integers, including zero.
  #
  R <- "full"

  #
  # Prior precision of zero-mean Gaussian prior for U and V of
  # the low-rank ARD prior.
  #
  lambda <- 0.1
  
  #
  # Whether to use the rotation explained in the ICML'11 paper.
  # Using the rotation is recommended, but for some data sets
  # equally good solutions can be obtained without rotation and
  # possibly faster.
  #  - TRUE|FALSE
  #
  rotate <- TRUE
  
  #
  # Parameters for controlling how the rotation is solved
  #  - opt.method chooses the optimization method and
  #    takes values "BFGS" or "L-BFGS". The former
  #    is typically faster but takes more memory, so the latter
  #    is the default choice. For small K may use BFGS instead.
  #  - opt.iter is the maximum number of iterations
  #  - lbfgs.factr is convergence criterion for L-BFGS; smaller
  #    values increase the accuracy (10^7 or 10^10 could be tried
  #    to speed things up)
  #  - bfgs.crit is convergence criterion for BFGS; smaller
  #    values increase the accuracy (10^-7 or 10^-3 could also be used)
  #
  opt.method <- "L-BFGS"
  opt.iter <- 10^5
  lbfgs.factr <- 10^10
  bfgs.crit <- 10^-5
  
  #
  # Initial value for the noise precisions. Should be large enough
  # so that the real structure is modeled with components
  # instead of the noise parameters (see Luttinen&Ilin, 2010)
  #  Values: Positive numbers, but generally should use values well
  #          above 1
  #
  init.tau <- 10^3
  
  #
  # Parameters for controlling when the algorithm stops.
  # It stops when the relative difference in the lower bound
  # falls below iter.crit or iter.max iterations have been performed.
  #
  iter.crit <- 10^-6
  iter.max <- 10^5
  
  #
  # Additive noise level for latent variables. The latent variables
  # of inactive components (those with very large alpha) occasionally
  # show some structure in the mean values, even though the distribution
  # matches very accurately the prior N(0,I). This structure disappears
  # is a tiny amount of random noise is added on top of the
  # mean estimates. Setting the value to 0 will make the predictions
  # deterministic
  #
  addednoise <- 1e-5
  
  #
  # Hyperparameters
  # - alpha_0, beta_0 for the ARD precisions
  # - alpha_0t, beta_0t for the residual noise predicions
  #
  prior.alpha_0 <- prior.beta_0 <- 1e-14
  prior.alpha_0t <- prior.beta_0t <- 1e-14
  
  # Two performace enhancing settings
  #   dropK=TRUE  : matrix dimensions are reduced when components are shut off
  #   low.mem=TRUE: list covW is not stored (may cause minor numerical differences)
  dropK <- TRUE
  low.mem <- FALSE
  
  #
  # Verbosity level
  #  0: Nothing
  #  1: Final cost function value for each run of GFAexperiment()
  #  2: Cost function values for each iteration
  #
  verbose <- 2
  
  return(list(R=R,rotate=rotate, init.tau=init.tau, iter.crit=iter.crit,
              iter.max=iter.max, opt.method=opt.method,
              lbfgs.factr=lbfgs.factr, bfgs.crit=bfgs.crit, opt.iter=opt.iter,
              addednoise=addednoise,
              prior.alpha_0=prior.alpha_0,prior.beta_0=prior.beta_0,
              prior.alpha_0t=prior.alpha_0t,prior.beta_0t=prior.beta_0t,
              verbose=verbose,lambda=lambda,dropK=dropK,low.mem=low.mem))
  
}
