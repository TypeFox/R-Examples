##########################################################################
## sample from a K-dimensional two-parameter item response model with
## logit link that has been rescaled so that the inverse link is:
##
## \delta0 + (1 - \delta0 - \delta1)*\Phi(.)
##
## where \delta0 \in (0, k0) and \delta1 \in (0, k1)
##
## priors for deltas are rescaled beta with parameters c0, d0, and c1, d1
##
##
## datamatrix is assumed to be nsubjects by nitems
##
## Andrew D. Martin
## Washington University
##
## Kevin M. Quinn
## Harvard University
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Feb. 17, 2005
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################



"MCMCirtKdRob" <-
  function(datamatrix, dimensions, item.constraints=list(),
           ability.constraints=list(),
           burnin = 500, mcmc = 5000, thin=1, interval.method="step",
           theta.w=0.5, theta.mp=4, alphabeta.w=1.0, alphabeta.mp=4,
           delta0.w=NA, delta0.mp=3, delta1.w=NA, delta1.mp=3,
           verbose = FALSE, seed = NA, theta.start=NA,
           alphabeta.start = NA, delta0.start = NA,
           delta1.start = NA, b0=0, B0=0,
           k0=.1, k1=.1, c0=1, d0=1, c1=1, d1=1,
           store.item=TRUE, store.ability=FALSE,
           drop.constant.items=TRUE, ... ) {
    
    ## set X up       
    if (drop.constant.items==TRUE){
      x.col.var <- apply(datamatrix, 2, var, na.rm=TRUE)
      keep.inds <- x.col.var>0
      keep.inds[is.na(keep.inds)] <- FALSE      
      datamatrix <- datamatrix[,keep.inds]
    }
    X <- as.data.frame(datamatrix)
    xvars <- dimnames(X)[[2]]
    xobs <- dimnames(X)[[1]]
    N <- nrow(X)    # number of subjects
    K <- ncol(X)    # number of items
    for (i in 1:K){
      X[is.na(X[,i]), i] <- -999
    }
    if(sum(datamatrix==1 | datamatrix==0 | is.na(datamatrix)) != (N*K)) {
      cat("Error: Data matrix contains elements other than 0, 1 or NA.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)
    }
    X <- as.matrix(X)

    ## take care of the case where X has no row names
    if (is.null(xobs)){
      xobs <- 1:N
    }
    
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)

    ## check slice sampling parameters
    if (!(interval.method %in% c("step", "doubling"))){
      cat("Error: interval.method not equal to 'step' or 'doubling'.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)      
    }
    method.step <- 0
    if (interval.method == "step"){
      method.step <- 1
    }
    if (theta.w <= 0 ){
      cat("Error: theta.w not > 0.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)      
    }
    if (theta.mp < 1 ){
      cat("Error: theta.mp not >= 1.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)      
    }
    if (alphabeta.w <= 0 ){
      cat("Error: alphabeta.w not > 0.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)      
    }
    if (alphabeta.mp < 1 ){
      cat("Error: alphabeta.mp not >= 1.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)      
    }
    if (is.na(delta0.w)){
      delta0.w <- 0.25*k0
    }
    if (delta0.w <= 0 ){
      cat("Error: delta0.w not > 0.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)      
    }
    if (delta0.mp < 1 ){
      cat("Error: delta0.mp not >= 1.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)      
    }
    if (is.na(delta1.w)){
      delta1.w <- 0.25*k1
    }
    if (delta1.w <= 0 ){
      cat("Error: delta1.w not > 0.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)      
    }
    if (delta1.mp < 1 ){
      cat("Error: delta1.mp not >= 1.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)      
    }
    
    ## error check the prior parameters for delta
    if (k0 < 0 | k0 > 0.5){
      cat("Error: k0 not in (0, 0.5).\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)      
    }
    if (k1 < 0 | k1 > 0.5){
      cat("Error: k1 not in (0, 0.5).\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)      
    }
    if (c0 < 0){
      cat("Error: c0 < 0.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)      
    }
    if (c1 < 0){
      cat("Error: c1 < 0.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)      
    }
    if (d0 < 0){
      cat("Error: d0 < 0.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)      
    }
    if (d1 < 0){
      cat("Error: d1 < 0.\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)      
    }
    
    ## setup constraints on Lambda = (alpha, beta)
    holder <- build.factor.constraints(item.constraints, X, K, dimensions+1)
    Lambda.eq.constraints <- holder[[1]]
    Lambda.ineq.constraints <- holder[[2]]
    X.names <- holder[[3]]

    ## setup constraints on theta
    holder <- build.factor.constraints(ability.constraints, t(X),
                                       N, dimensions)
    theta.eq.constraints <- holder[[1]]
    theta.ineq.constraints <- holder[[2]]

    ## setup prior on Lambda
    holder <- form.factload.norm.prior(b0, B0, K, dimensions+1, X.names)
    Lambda.prior.mean <- holder[[1]]
    Lambda.prior.prec <- holder[[2]]

    # seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## Starting values for delta0 and delta1
    if (is.na(delta0.start)){
      delta0.start <- 0.5 * k0;
    }
    if (is.na(delta1.start)){
      delta1.start <- 0.5 * k1;
    }
    if (delta0.start < 0 | delta0.start > k0){
      cat("Error: delta0 not in (0, k0).\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)      
    }
    if (delta1.start < 0 | delta1.start > k1){
      cat("Error: delta1 not in (0, k1).\n")
      stop("Please check data and try MCMCirtKdRob() again.\n",
           call.=FALSE)      
    }

    
    
    ## Starting values for Lambda
    Lambda <- matrix(0, K, dimensions+1)
    if (is.na(alphabeta.start)){# sets Lambda to equality constraints & 0s
      for (i in 1:K){
        for (j in 1:(dimensions+1)){
          if (Lambda.eq.constraints[i,j]==-999){
            if(Lambda.ineq.constraints[i,j]==0){
              if (j==1){
                  probit.out <- glm(as.factor(X[X[,i]!=-999,i])~1,
                                    family=binomial(link="logit"))
                  probit.beta <- coef(probit.out)
                  Lambda[i,j] <- -1 * probit.beta[1]
              }
            }
            if(Lambda.ineq.constraints[i,j]>0){
              Lambda[i,j] <- 1.0
            }
            if(Lambda.ineq.constraints[i,j]<0){
              Lambda[i,j] <- -1.0
            }          
          }
          else Lambda[i,j] <- Lambda.eq.constraints[i,j]
        }
      }    
    }
    else if (is.matrix(alphabeta.start)){
      if (nrow(alphabeta.start)==K && ncol(alphabeta.start)==(dimensions+1))
        Lambda  <- alphabeta.start
      else {
        cat("Starting values not of correct size for model specification.\n")
        stop("Please respecify and call MCMCirtKdRob() again\n", call.=FALSE)
      }
    }
    else if (length(alphabeta.start)==1 && is.numeric(alphabeta.start)){
      Lambda <- matrix(alphabeta.start, K, dimensions+1)
      for (i in 1:K){
        for (j in 1:(dimensions+1)){
          if (Lambda.eq.constraints[i,j] != -999)
            Lambda[i,j] <- Lambda.eq.constraints[i,j]
        }
      }    
    }
    else {
      cat("Starting values for alpha & beta neither NA, matrix, nor scalar.\n")
      stop("Please respecify and call MCMCirtKdRob() again\n", call.=FALSE)
    }
    for (i in 1:K){
      lam.sqdist <- sum(Lambda[i,]^2)
      while (lam.sqdist > 100){
        Lambda[i,] <- Lambda[i,] * 0.95
        lam.sqdist <- sum(Lambda[i,]^2)        
      }
    }
    
    
    
    ## Starting values for theta
    if (is.na(theta.start)){
      theta <- matrix(0, N, dimensions)
    }
    else if(is.null(dim(theta.start)) & is.numeric(theta.start)){
      theta <- matrix(theta.start, N, dimensions)
    }
    else if(nrow(theta.start)==N & ncol(theta.start)==dimensions){
      theta <- theta.start
    }
    else{
      cat("Starting values for theta neither NA, matrix, nor scalar.\n")
      stop("Please respecify and call MCMCirtKdRob() again\n", call.=FALSE)     
    }
    for (i in 1:N){
      for (j in 1:dimensions){
        if (theta.eq.constraints[i,j]==-999){
          if(theta.ineq.constraints[i,j]>0){
            theta[i,j] <- 0.5
          }
          if(theta.ineq.constraints[i,j]<0){
            theta[i,j] <- -0.5
          }          
        }
        else theta[i,j] <- theta.eq.constraints[i,j]
      }
    }
      
    


    
    ## define holder for posterior sample
    if (store.ability == FALSE && store.item == FALSE){
      cat("You need to store either the ability or item parameters.\n")
      stop("Please respecify and call MCMCirtKdRob() again\n", call.=FALSE)  
    }
    else if (store.ability == TRUE && store.item == FALSE){
      sample <- matrix(data=0, mcmc/thin, (dimensions+1)*N+2)
    }
    else if(store.ability == FALSE && store.item == TRUE) {
      sample <- matrix(data=0, mcmc/thin, K*(dimensions+1)+2)
    }
    else { # store.ability==TRUE && store.item==TRUE
      sample <- matrix(data=0, mcmc/thin, K*(dimensions+1)+(dimensions+1)*N+2)
    }

        
    ## Call the C++ code to do the real work
    posterior <- .C("irtKdRobpost",
                    samdata = as.double(sample),
                    samrow = as.integer(nrow(sample)),
                    samcol = as.integer(ncol(sample)),
                    X = as.integer(X),
                    Xrow = as.integer(nrow(X)),
                    Xcol = as.integer(ncol(X)),
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                    verbose = as.integer(verbose),
                    method.step = as.integer(method.step),
                    theta.w = as.double(theta.w),
                    theta.mp = as.integer(theta.mp),
                    ab.w = as.double(alphabeta.w),
                    ab.mp = as.integer(alphabeta.mp),
                    delta0.w = as.double(delta0.w),
                    delta0.mp = as.integer(delta0.mp),
                    delta1.w = as.double(delta1.w),
                    delta1.mp = as.integer(delta1.mp),                    
                    delta0 = as.double(delta0.start),
                    delta1 = as.double(delta1.start),
                    Lambda = as.double(Lambda),
                    Lambdarow = as.integer(nrow(Lambda)),
                    Lambdacol = as.integer(ncol(Lambda)),
                    theta = as.double(theta),
                    thetarow = as.integer(nrow(theta)),
                    thetacol = as.integer(ncol(theta)),                    
                    Lameq = as.double(Lambda.eq.constraints),
                    Lameqrow = as.integer(nrow(Lambda.eq.constraints)),
                    Lameqcol = as.integer(ncol(Lambda.ineq.constraints)),
                    Lamineq = as.double(Lambda.ineq.constraints),
                    Lamineqrow = as.integer(nrow(Lambda.ineq.constraints)),
                    Lamineqcol = as.integer(ncol(Lambda.ineq.constraints)),
                    theteq = as.double(theta.eq.constraints),
                    theteqrow = as.integer(nrow(theta.eq.constraints)),
                    theteqcol = as.integer(ncol(theta.ineq.constraints)),
                    thetineq = as.double(theta.ineq.constraints),
                    thetineqrow = as.integer(nrow(theta.ineq.constraints)),
                    thetineqcol = as.integer(ncol(theta.ineq.constraints)),
                    Lampmean = as.double(Lambda.prior.mean),
                    Lampmeanrow = as.integer(nrow(Lambda.prior.mean)),
                    Lampmeancol = as.integer(ncol(Lambda.prior.prec)),
                    Lampprec = as.double(Lambda.prior.prec),
                    Lampprecrow = as.integer(nrow(Lambda.prior.prec)),
                    Lamppreccol = as.integer(ncol(Lambda.prior.prec)),
                    k0 = as.double(k0),
                    k1 = as.double(k1),
                    c0 = as.double(c0),
                    c1 = as.double(c1),
                    d0 = as.double(d0),
                    d1 = as.double(d1),
                    storeitem = as.integer(store.item),
                    storesability = as.integer(store.ability),
                    PACKAGE="MCMCpack"
                    )


    
    ## put together matrix and build MCMC object to return
    sample <- matrix(posterior$samdata, posterior$samrow, posterior$samcol,
                     byrow=FALSE)
    output <- mcmc(data=sample,start=1, end=mcmc, thin=thin)
    
    par.names <- NULL
    if (store.item==TRUE){
      alpha.hold <- paste("alpha", X.names, sep=".")
      beta.hold <- paste("beta", X.names, sep = ".")
      beta.hold <- rep(beta.hold, dimensions, each=dimensions)
      beta.hold <- paste(beta.hold, 1:dimensions, sep=".")
      
      Lambda.names <- t(cbind(matrix(alpha.hold, K, 1), 
                              matrix(beta.hold,K,dimensions,byrow=TRUE)))  
      dim(Lambda.names) <- NULL
      par.names <- c(par.names, Lambda.names)
    }
    
    if (store.ability==TRUE){
      phi.names <- paste(paste("theta",
                               rep(xobs, each=(dimensions+1)), sep="."),
                         rep(0:dimensions,(dimensions+1)), sep=".")
      par.names <- c(par.names, phi.names)            
    }
    
    par.names <- c("delta0", "delta1", par.names)

    varnames(output) <- par.names

    ## get rid of columns for constrained parameters
    output.df <- as.data.frame(as.matrix(output))
    output.var <- diag(var(output.df))
    output.df <- output.df[,output.var != 0]
    output <- mcmc(as.matrix(output.df), start=1, end=mcmc, thin=thin)
    
    ## add constraint info so this isn't lost
    attr(output, "constraints") <- item.constraints
    attr(output, "n.items") <- K
    attr(output, "n.dimensions") <- dimensions
    attr(output,"title") <-
      "MCMCpack Robust K-Dimensional Item Response Theory Model Posterior Sample" 
    
return(output)

    
  }
