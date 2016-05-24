########## hidden functions to help in model implementation ##########

## NOTE: these are not exported to the user and should always be
##       used in model functions.  As such, fixing problems here
##       fixes them in all functions simultaneously.
##
## updated by ADM 7/22/04
## re-organized (alphabetical) by ADM 7/28/04
## added a number of functions for teaching models by ADM 1/25/2006 

## create an agreement score matrix from a vote matrix
## subjects initially on rows and items on cols of X
## note: treats missing votes as category for agreement / might be
##       more principled to treat them in another fashion
"agree.mat" <- function(X){
  X <- t(X) # put subjects on columns
  n <- ncol(X)
  X[is.na(X)] <- -999
  A <- matrix(NA, n, n)
  for (i in 1:n){
    A[i,] <- apply(X[,i] == X, 2, sum)
  }
  return(A/nrow(X))
}

## create constraints for measurement models
"build.factor.constraints" <-
  function(lambda.constraints, X, K, factors){

    ## build initial constraint matrices and assign var names
    Lambda.eq.constraints <- matrix(NA, K, factors)
    Lambda.ineq.constraints <- matrix(0, K, factors)    
    if (is.null(colnames(X))){
      X.names <- paste("V", 1:ncol(X), sep="")
    }
    else {
      X.names <- colnames(X)
    }    
    rownames(Lambda.eq.constraints) <- X.names
    rownames(Lambda.ineq.constraints) <- X.names

    ## setup the equality and inequality contraints on Lambda
    if (length(lambda.constraints) != 0){
      constraint.names <- names(lambda.constraints)  
      for (i in 1:length(constraint.names)){
        name.i <- constraint.names[i]
        lambda.constraints.i <- lambda.constraints[[i]]
        col.index <- lambda.constraints.i[[1]]
        replace.element <- lambda.constraints.i[[2]]
        if (is.numeric(replace.element)){
          Lambda.eq.constraints[rownames(Lambda.eq.constraints)==name.i,
                                col.index] <- replace.element
        }
        if (replace.element=="+"){
          Lambda.ineq.constraints[rownames(Lambda.ineq.constraints)==name.i,
                                  col.index] <- 1 
        }
        if (replace.element=="-"){
          Lambda.ineq.constraints[rownames(Lambda.ineq.constraints)==name.i,
                                  col.index] <- -1
        }
      }
    }
    
    testmat <- Lambda.ineq.constraints * Lambda.eq.constraints
    
    if (min(is.na(testmat))==0){
      if ( min(testmat[!is.na(testmat)]) < 0){
        cat("Constraints on factor loadings are logically inconsistent.\n")
        stop("Please respecify and call ", calling.function(), " again.\n")
      }
    }
    Lambda.eq.constraints[is.na(Lambda.eq.constraints)] <- -999
    
    return( list(Lambda.eq.constraints, Lambda.ineq.constraints, X.names))    
  }

# return name of the calling function
"calling.function" <-
   function(parentheses=TRUE) {
     calling.function <- strsplit(toString(sys.call(which=-3)),",")[[1]][1]
     if (parentheses){
       calling.function <- paste(calling.function, "()", sep="")
     }
     return(calling.function)
   }

# check inverse Gamma prior
"check.ig.prior" <-
   function(c0, d0) {
     
    if(c0 <= 0) {
      cat("Error: IG(c0/2,d0/2) prior c0 less than or equal to zero.\n")
      stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
    }
    if(d0 <= 0) {
      cat("Error: IG(c0/2,d0/2) prior d0 less than or equal to zero.\n")
       stop("Please respecify and call ", calling.function(), " again.\n",
        call.=FALSE)    
    }
    return(0)
  }

# check beta prior
"check.beta.prior" <-
   function(alpha, beta) {
     
    if(alpha <= 0) {
      cat("Error: Beta(alpha,beta) prior alpha less than or equal to zero.\n")
      stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
    }
    if(beta <= 0) {
      cat("Error: Beta(alpha,beta) prior beta less than or equal to zero.\n")
       stop("Please respecify and call ", calling.function(), " again.\n",
        call.=FALSE)    
    }
    return(0)
  }

# check Gamma prior
# ADM 1/25/2006
"check.gamma.prior" <-
   function(alpha, beta) {
     
    if(alpha <= 0) {
      cat("Error: Gamma(alpha,beta) prior alpha less than or equal to zero.\n")
      stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
    }
    if(alpha <= 0) {
      cat("Error: Gamma(alpha,beta) prior beta less than or equal to zero.\n")
       stop("Please respecify and call ", calling.function(), " again.\n",
        call.=FALSE)    
    }
    return(0)
  }   

# check Normal prior
# ADM 1/26/2006
"check.normal.prior" <-
   function(mu, sigma2) {
   
     if(sigma2 <= 0) {
       cat("Error: Normal(mu0,tau20) prior sigma2 less than or equal to zero.\n")
       stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)    
     }
   }

# check mc parameter
# ADM 1/25/2006
"check.mc.parameter" <-
  function(mc) {
  
    if(mc < 0) {
      cat("Error: Monte Carlo iterations negative.\n")
      stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE) 
    }   
    return(0)
  }

# check mcmc parameters
"check.mcmc.parameters" <-
  function(burnin, mcmc, thin) {
  
    if(mcmc %% thin != 0) {
      cat("Error: MCMC iterations not evenly divisible by thinning interval.\n")
      stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
    }
    if(mcmc < 0) {
      cat("Error: MCMC iterations negative.\n")
      stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE) 
    }
    if(burnin < 0) {
      cat("Error: Burnin iterations negative.\n")
      stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
    }
    if(thin < 0) {
      cat("Error: Thinning interval negative.\n")
      stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
    }    
    return(0)
  }

# check to see if an offset is passed
"check.offset" <-
   function(args) {
      if(sum(names(args)=="offset")==1) {
         cat("Error: Offsets are currently not supported in MCMCpack.\n")
         stop("Please respecify and call ", calling.function(), " again.\n",
            call.=FALSE)
      }
   return(0)
   }

# put together starting values for coefficients
# NOTE: This can be used for any GLM model by passing the right family
#       or for another model by passing default starting values to
#       the function   
"coef.start" <-
   function(beta.start, K, formula, family, data=NULL, defaults=NA) {
          
     if (is.na(beta.start)[1] & is.na(defaults)[1]){ # use GLM estimates
       beta.start <- matrix(coef(glm(formula, family=family, data=data)), K, 1)
     }
     else if(is.na(beta.start)[1] & !is.na(defaults)[1]){ # use passed values
        beta.start <- matrix(defaults,K,1)     
     }
     else if(is.null(dim(beta.start))) {
       beta.start <- beta.start * matrix(1,K,1)  
       }
     else if(!all(dim(beta.start) == c(K,1))) {
       cat("Error: Starting value for beta not conformable.\n")
       stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
       }
     return(beta.start)
   }

## generate starting values for a factor loading matrix
"factload.start" <-
  function(lambda.start, K, factors, Lambda.eq.constraints,
           Lambda.ineq.constraints){

    Lambda <- matrix(0, K, factors)
    if (any(is.na(lambda.start))){# sets Lambda to equality constraints & 0s
      for (i in 1:K){
        for (j in 1:factors){
          if (Lambda.eq.constraints[i,j]==-999){
            if(Lambda.ineq.constraints[i,j]==0){
              Lambda[i,j] <- 0
            }
            if(Lambda.ineq.constraints[i,j]>0){
              Lambda[i,j] <- .5
            }
            if(Lambda.ineq.constraints[i,j]<0){
              Lambda[i,j] <- -.5
            }          
          }
          else Lambda[i,j] <- Lambda.eq.constraints[i,j]
        }
      }    
    }
    else if (is.matrix(lambda.start)){
      if (nrow(lambda.start)==K && ncol(lambda.start)==factors)
        Lambda  <- lambda.start
      else {
        cat("lambda.start not of correct size for model specification.\n")
        stop("Please respecify and call ", calling.function(), " again.\n")
      }
    }
    else if (length(lambda.start)==1 && is.numeric(lambda.start)){
      Lambda <- matrix(lambda.start, K, factors)
      for (i in 1:K){
        for (j in 1:factors){
          if (Lambda.eq.constraints[i,j] != -999)
            Lambda[i,j] <- Lambda.eq.constraints[i,j]
        }
      }    
    }
    else {
      cat("lambda.start neither NA, matrix, nor scalar.\n")
      stop("Please respecify and call ", calling.function, " again.\n")
    }

    return(Lambda)
  }




## based on code originally written by Keith Poole
## takes a subject by subject agreement score matrix as input
"factor.score.eigen.start" <- function(A, factors){

  A <- (1 - A)^2
  AA <- A
  arow <- matrix(NA, nrow(A), 1)
  acol <- matrix(NA, ncol(A), 1)
  
  for (i in 1:nrow(A)){
    arow[i] <- mean(A[i,])
  }
  for (i in 1:ncol(A)){
    acol[i] <- mean(A[,i])
  }
  
  matrixmean <- mean(acol)
  
  for (i in 1:nrow(A)){
    for (j in 1:ncol(A)){
      AA[i,j] <- (A[i,j]-arow[i]-acol[j]+matrixmean)/(-2)
    }
  }
  
  ev <- eigen(AA)
  scores <- matrix(NA, nrow(A), factors)
  for (i in 1:factors){
    scores[,i] <- ev$vec[,i]*sqrt(ev$val[i])
    scores[,i] <- (scores[,i] - mean(scores[,i]))/sd(scores[,i])
  }
  return(scores)
}


## check starting values of factor scores or ability parameters
## subjects on rows of X
"factor.score.start.check" <- function(theta.start, X, prior.mean,
                                       prior.prec, eq.constraints,
                                       ineq.constraints, factors){

  N <- nrow(X)
  
  ## set value of theta.start
  if (max(is.na(theta.start))==1) { 
    theta.start <- factor.score.eigen.start(agree.mat(X), 1)
    for (i in 1:factors){
      theta.start[,i] <- prior.mean[i] + theta.start[,i] *
        sqrt(1/prior.prec[i,i])
      
      # make sure these are consistent with hard and soft constraints  
      for (j in 1:nrow(theta.start)){
        if (eq.constraints[j,i] != -999){
          if (theta.start[j,i] * eq.constraints[j,i] < 0){
            theta.start[,i] <- -1*theta.start[,i]
          }
        }
        if (theta.start[j,i] * ineq.constraints[j,i] < 0){
          theta.start[,i] <- -1*theta.start[,i]
        }        
      }
      theta.start[eq.constraints[,i]!=-999,i] <-
         eq.constraints[eq.constraints[,i]!=-999,i]
      theta.start[ineq.constraints[,i]!=0,i] <-
         abs(theta.start[ineq.constraints[,i]!=0,i]) * 
         ineq.constraints[ineq.constraints[,i]!=0,i]
    }
  }
  else if(is.numeric(theta.start) && is.null(dim(theta.start))) {
    theta.start <- theta.start * matrix(1, N, 1)  
  }
  else if((dim(theta.start)[1] != N) ||
           (dim(theta.start)[2] != factors)) {
    cat("Starting value for theta not appropriately sized.\n")
    stop("Please respecify and call", calling.function(), " again.\n",
         call.=FALSE)
  }
  else {
    cat("Inappropriate value of theta.start passed.\n")
    stop("Please respecify and call", calling.function(), " again.\n",
         call.=FALSE)      
  }

  ## check value of theta.start
  prev.bind.constraints <- rep(0, factors)
  for (i in 1:N){
    for (j in 1:factors){
      if (eq.constraints[i,j]==-999){
        if(ineq.constraints[i,j]>0 && theta.start[i,j] < 0){
          if (prev.bind.constraints[j]==0){
            theta.start[,j] <- -1*theta.start[,j]
          }
          else {
            cat("Parameter constraints logically inconsistent.\n")
            stop("Please respecify and call ", calling.function(), " again.",
                 call.=FALSE)                
          }
          prev.bind.constraints[j] <- prev.bind.constraints[j] + 1
        }
        if(ineq.constraints[i,j]<0 && theta.start[i,j] > 0){
          if (prev.bind.constraints[j]==0){
            theta.start[,j] <- -1*theta.start[,j]
          }
          else {
            cat("Parameter constraints logically inconsistent.\n")
            stop("Please respecify and call ", calling.function(), " again.",
                 call.=FALSE)                
          }
          prev.bind.constraints[j] <- prev.bind.constraints[j] + 1
        }
      }
      else {
        if ((theta.start[i,j] * eq.constraints[i,j]) > 0){
          theta.start[i,j] <- eq.constraints[i,j]
        }
        else {
          if (prev.bind.constraints[j]==0){
            theta.start[,j] <- -1*theta.start[,j]
            theta.start[i,j] <- eq.constraints[i,j]
          }
          else {
            cat("Parameter constraints logically inconsistent.\n")
            stop("Please respecify and call ", calling.function(), " again.",
                 call.=FALSE)                
          }
          prev.bind.constraints[j] <- prev.bind.constraints[j] + 1
        }
      }      
    }
  }    
  return(theta.start)
}

## get starting values for factor uniqueness matrix (Psi)
"factuniqueness.start" <-
  function(psi.start, X){
    
    K <- ncol(X)
    if (any(is.na(psi.start))){
      Psi <- 0.5 * diag(diag(var(X)))
    }
    else if (is.double(psi.start) &&
             (length(psi.start==1) || length(psi.start==K))){
      Psi <- diag(K) * psi.start
    }
    else {
      cat("psi.start neither NA, double. nor appropriately sized matrix.\n")
      stop("Please respecify and call ", calling.function, " again.\n")
    }
    if (nrow(Psi) != K || ncol(Psi) != K){
      cat("Psi starting value not K by K matrix.\n")
      stop("Please respecify and call ", calling.function, " again.\n")    
    }

    return(Psi)
  }


## form the ind. normal prior for a factor loading matrix
"form.factload.norm.prior" <-
  function(l0, L0, K, factors, X.names){

    ## prior means
    if (is.matrix(l0)){ # matrix input for l0
      if (nrow(l0)==K && ncol(l0)==factors){
        Lambda.prior.mean <- l0
        rownames(Lambda.prior.mean) <- X.names
      }
      else {
        cat("l0 not of correct size for model specification.\n")
        stop("Please respecify and call ", calling.function(), " again.\n")
      }
    }
    else if (is.list(l0)){ # list input for l0
      Lambda.prior.mean <- matrix(0, K, factors)
      rownames(Lambda.prior.mean) <- X.names
      l0.names <- names(l0)
      for (i in 1:length(l0.names)){
        name.i <- l0.names[i]
        l0.i <- l0[[i]]
        col.index <- l0.i[[1]]
        replace.element <- l0.i[[2]]
        if (is.numeric(replace.element)){
          Lambda.prior.mean[rownames(Lambda.prior.mean)==name.i,
                            col.index] <- replace.element
        }   
      }
    }
    else if (length(l0)==1 && is.numeric(l0)){ # scalar input for l0
      Lambda.prior.mean <- matrix(l0, K, factors)
      rownames(Lambda.prior.mean) <- X.names
    }
    else {
      cat("l0 neither matrix, list, nor scalar.\n")
      stop("Please respecify and call ", calling.function(), " again.\n")
    }
    
    ## prior precisions
    if (is.matrix(L0)){ # matrix input for L0
      if (nrow(L0)==K && ncol(L0)==factors){
        Lambda.prior.prec <- L0
        rownames(Lambda.prior.prec) <- X.names
      }
      else {
        cat("L0 not of correct size for model specification.\n")
        stop("Please respecify and call ", calling.function(), " again.\n")
      }
    }
    else if (is.list(L0)){ # list input for L0
      Lambda.prior.prec <- matrix(0, K, factors)  
      rownames(Lambda.prior.prec) <- X.names
      L0.names <- names(L0)
      for (i in 1:length(L0.names)){
        name.i <- L0.names[i]
        L0.i <- L0[[i]]
        col.index <- L0.i[[1]]
        replace.element <- L0.i[[2]]
        if (is.numeric(replace.element)){
          Lambda.prior.prec[rownames(Lambda.prior.prec)==name.i,
                            col.index] <- replace.element
        }   
      }
    }
    else if (length(L0)==1 && is.numeric(L0)){ # scalar input for L0
      Lambda.prior.prec <- matrix(L0, K, factors)
      rownames(Lambda.prior.prec) <- X.names
    }
    else {
      cat("L0 neither matrix, list, nor scalar.\n")
      stop("Please respecify and call ", calling.function(), " again.\n")
    }
    if (min(Lambda.prior.prec) < 0) {
      cat("L0 contains negative elements.\n")
      stop("Please respecify and call ", calling.function(), " again.\n")
    } 


    return( list(Lambda.prior.mean, Lambda.prior.prec))    
  }

## form ind. inv. gamma prior for a diagonal var. cov. matrix
"form.ig.diagmat.prior" <-
  function(a0, b0, K){

    ## setup prior for diag(Psi)
    if (length(a0)==1 && is.double(a0))
      a0 <- matrix(a0, K, 1)
    else if (length(a0) == K && is.double(a0))
      a0 <- matrix(a0, K, 1)
    else {
      cat("a0 not properly specified.\n")
      stop("Please respecify and call ", calling.function, " again.\n")
    }
    if (length(b0)==1 && is.double(b0))
      b0 <- matrix(b0, K, 1)
    else if (length(b0) == K && is.double(b0))
      b0 <- matrix(b0, K, 1)
    else {
      cat("b0 not properly specified.\n")
      stop("Please respecify and call ", calling.function(), " again.\n")
    }
    
    ## prior for Psi error checking
    if(min(a0) <= 0) {
      cat("IG(a0/2,b0/2) prior parameter a0 less than or equal to zero.\n")
      stop("Please respecify and call ", calling.function, " again.\n")
    }
    if(min(b0) <= 0) {
      cat("IG(a0/2,b0/2) prior parameter b0 less than or equal to zero.\n")
      stop("Please respecify and call ", calling.function(), " again.\n")      
    }  

    return(list(a0, b0) )
  }

# pull together the posterior density sample
"form.mcmc.object" <-
  function(posterior.object, names, title, ...) {
    holder <- matrix(posterior.object$sampledata,
                     posterior.object$samplerow,
                     posterior.object$samplecol,
                     byrow=FALSE)
      
    output <- mcmc(data=holder, start=(posterior.object$burnin+1),
                   end=(posterior.object$burnin+posterior.object$mcmc),
                   thin=posterior.object$thin)
    varnames(output) <- as.list(names)
    attr(output,"title") <- title
    
    attribs <- list(...)
    K <- length(attribs)
    attrib.names <- names(attribs)

    if (K>0){
      for (i in 1:K){
        attr(output, attrib.names[i]) <- attribs[[i]]
      }
    }
    
    return(output)  
  }

# form multivariate Normal prior
"form.mvn.prior" <-
   function(b0, B0, K) {
  
     # prior mean
     if(is.null(dim(b0))) {
       b0 <- b0 * matrix(1,K,1)  
     } 
     if((dim(b0)[1] != K) || (dim(b0)[2] != 1)) {
       cat("Error: N(b0,B0^-1) prior b0 not conformable.\n")
       stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
     }
     
     # prior precision
     if(is.null(dim(B0))) {
       if (length(B0) > K){
         stop("B0 was passed as a vector longer than K.\nB0 must be either a scalar or a matrix.\nPlease respecify and call ", calling.function(), " again.\n", call.=FALSE)
       }
       B0 <- B0 * diag(K)    
     }
     if((dim(B0)[1] != K) || (dim(B0)[2] != K)) {
       cat("Error: N(b0,B0^-1) prior B0 not conformable.\n")
       stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
     }

     ## check B0 for symmetry
     symproblem <- FALSE
     for (i in 1:K){
       for (j in i:K){
         if (B0[i,j] != B0[j,i]){
           symproblem <- TRUE
         }
       }
     }
     if (symproblem){
       cat("B0 is not symmetric.\n")
       stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)       
     }
     
     return(list(b0,B0))
   }

# parse the passed seeds
# 1] if a scalar is passed, it is used by Mersennse twister
# 2] if a list of length two is passed, a parallel-friendly stream is
#    created using L'Ecuyer
"form.seeds" <-
   function(seed) {
      if(length(seed)==1) {
         if(is.na(seed)) seed <- 12345
         seed <- as.integer(seed)
         if(seed < 0) {
            cat("Error: Mersenne seed negative.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)                       
         }
         seeds <- list(0, rep(seed,6), 0)
      }
      if(length(seed)==2) {
         if(!is.list(seed)) {
            cat("Error: List must be passed to use L'Ecuyer.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)          
         }
         lec.seed <- seed[[1]]
         lec.substream <- as.integer(seed[[2]])
         if(is.na(lec.seed[1])) lec.seed <- rep(12345, 6)
         if(length(lec.seed) != 6) {
            cat("Error: L'Ecuyer seed not of length six.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)          
         }
         if(!all(lec.seed >= 0))  {
             cat("Error: At least one L'Ecuyer seed negative.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)          
         }
         if( max(lec.seed[1:3]) >= 4294967087){
           cat("Error: At least one of first three L'Ecuyer seeds\n")
           cat("  greater than or equal to 4294967087\n")
           stop("Please respecify and call ", calling.function(), " again.",
                call.=FALSE)          
         }
         if( all(lec.seed[1:3] == 0 )){
           cat("Error: first three L'Ecuyer seeds == 0\n")
           stop("Please respecify and call ", calling.function(), " again.",
                call.=FALSE)          
         }
         if( max(lec.seed[4:6]) >= 4294944443){
           cat("Error: At least one of last three L'Ecuyer seeds\n")
           cat("  greater than or equal to 4294944443\n")
           stop("Please respecify and call ", calling.function(), " again.",
                call.=FALSE)          
         }         
         if( all(lec.seed[4:6] == 0 )){
           cat("Error: last three L'Ecuyer seeds == 0\n")
           stop("Please respecify and call ", calling.function(), " again.",
                call.=FALSE)          
         }
         if(lec.substream < 1) {
            cat("Error: L'Ecuyer substream number not positive.\n")
            stop("Please respecify and call ", calling.function(), " again.",
                 call.=FALSE)               
         }
         seeds <- list(1, lec.seed, lec.substream) 
      }
      if(length(seed)>2) {
            cat("Error: Seed passed as length greater than two.\n")
            stop("Please respecify and call ", calling.function(), " again.",
              call.=FALSE)        
      }
      return(seeds)
   }

# form Wishart prior
"form.wishart.prior" <-
    function(v, S, K) {
    
    # check to see if degrees of freedom produces proper prior
    if(v < K) {
      cat("Error: Wishart(v,S) prior v less than or equal to K.\n")
      stop("Please respecify and call ", calling.function(), " again.\n")
    } 
    
    # form the prior scale matrix
    if(is.null(dim(S))) {
      S <- S * diag(K)
    }
    if((dim(S)[1] != K) | (dim(S)[2] != K)) {
      cat("Error: Wishart(v,S) prior S not comformable [K times K].\n")
      stop("Please respecify and call ", calling.function(), " again.\n")
    }
           
    return(list(v,S))
}

# parse formula and return a list that contains the model response
# matrix as element one, and the model matrix as element two
"parse.formula" <- 
   function(formula, data=NULL, intercept=TRUE, justX=FALSE) {

    # extract Y, X, and variable names for model formula and frame
    mf <- match.call(expand.dots = FALSE)
    mf$intercept <- mf$justX <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    mt <- attr(mf, "terms")
    if (!intercept){
      attributes(mt)$intercept <- 0
    }

    # null model support
    X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
    X <- as.matrix(X)         # X matrix
    xvars <- dimnames(X)[[2]] # X variable names
    xobs  <- dimnames(X)[[1]] # X observation names
    if (justX){
      Y <- NULL
    }
    else {
      Y <- as.matrix(model.response(mf, "numeric")) # Y matrix
    }
    return(list(Y, X, xvars, xobs))
   }

# setup tuning constant for scalar parameter
"scalar.tune" <- function(mcmc.tune){
  if (max(is.na(mcmc.tune))){
    cat("Error: Scalar tuning parameter cannot contain NAs.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)    
  }
  if (length(mcmc.tune) != 1){
    cat("Error: Scalar tuning parameter does not have length = 1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (mcmc.tune <= 0) {
    cat("Error: Scalar tuning parameter not positive.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  } 
  return(mcmc.tune)
}
  
# put together starting values for sigma2
"sigma2.start" <-
   function(sigma2.start, formula, data) {
   
     if(is.na(sigma2.start)){ # use MLE
       lm.out <- lm(formula, data=data)
       sigma2.start <- var(residuals(lm.out))
     }   
     else if(sigma2.start <= 0) {
       cat("Error: Starting value for sigma2 negative.\n")
       stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
     }
     else if (length(sigma2.start) != 1){
       cat("Error: Starting value for sigma2 not a scalar.\n")
       stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
     }
     else if (!is.numeric(sigma2.start)){
       cat("Error: Starting value for sigma2 neither numeric nor NA.\n")
       stop("Please respecify and call ", calling.function(), " again.\n",
         call.=FALSE)
     }
     return(sigma2.start)   
   }




## setup diagonal tuning matrix for vector parameters
"vector.tune" <- function(mcmc.tune, K){
  if (max(is.na(mcmc.tune))){
    cat("Error: Vector tuning parameter cannot contain NAs.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)    
  }
  if (length(mcmc.tune) == 1){
    mcmc.tune <- rep(mcmc.tune, K)
  }
  if (length(mcmc.tune) != K){
    cat("Error: length(vector tuning parameter) != length(theta) or 1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(mcmc.tune <= 0) != 0) {
    cat("Error: Vector tuning parameter cannot contain negative values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (length(mcmc.tune)==1){
    return(matrix(mcmc.tune, 1, 1))
  }
  else{
    return(diag(as.double(mcmc.tune)))
  }
}

