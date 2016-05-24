###########################################################################
## sample from the posterior distribution of a factor analysis model
## model in R using linked C++ code in Scythe.
##
## The model is:
##
## x*_i = \Lambda \phi_i + \epsilon_i,   \epsilon_i \sim N(0, I)
##
## \lambda_{ij} \sim N(l0_{ij}, L0^{-1}_{ij})
## \phi_i \sim N(0,I)
##
## and x*_i is the latent variable formed from the observed ordinal
## variable in the usual (Albert and Chib, 1993) way.
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
## May 12, 2003
## Revised to accommodate new spec 7/13/2004
## 
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

"MCMCordfactanal" <-
  function(x, factors, lambda.constraints=list(),
           data=parent.frame(), burnin = 1000, mcmc = 20000,
           thin=1, tune=NA, verbose = 0, seed = NA,
           lambda.start = NA, l0=0, L0=0,
           store.lambda=TRUE, store.scores=FALSE,
           drop.constantvars=TRUE, ... ) {

    ## check for MCMCirtKd special case, this is used to tell the R
    ## and C++ code what to echo (1 if regular, 2 if MCMCirtKd)
    ## the test is based on the existence of model="MCMCirtKd"
    ## passed through ...
    args <- list(...)

    if (length(args$model) > 0){ # if model arg is passed
      if (args$model=="MCMCirtKd"){
        case.switch <- 2
        echo.name <- "MCMCirtKd"
      }
      ## could allow for other possibities here but not clear what they
      ## would be 
    }
    else { # if model arg not passed then assume MCMCordfactanal
      case.switch <- 1
      echo.name <- "MCMCordfactanal"
    }
    
 
    # extract X and variable names from the model formula and frame       
    if (is.matrix(x)){
      if (drop.constantvars==TRUE){
        x.col.var <- apply(x, 2, var, na.rm=TRUE)
        keep.inds <- x.col.var>0
        keep.inds[is.na(keep.inds)] <- FALSE      
        x <- x[,keep.inds]
      }
      X <- as.data.frame(x)
      xvars <- dimnames(X)[[2]]
      xobs <- dimnames(X)[[1]]
      N <- nrow(X)    # number of observations
      K <- ncol(X)    # number of manifest variables
      ncat <- matrix(NA, K, 1) # vector of number of categ. in each man. var. 
      for (i in 1:K){
        X[,i] <- factor(X[,i], ordered=TRUE)
        ncat[i] <- nlevels(X[,i])
        X[,i] <- as.integer(X[,i])
        X[is.na(X[,i]), i] <- -999
      }
      X <- as.matrix(X)
    }
    else {
      call <- match.call()
      mt <- terms(x, data=data)
      if (attr(mt, "response") > 0) 
        stop("Response not allowed in formula in ", echo.name, "().\n")
      if(missing(data)) data <- sys.frame(sys.parent())
      mf <- match.call(expand.dots = FALSE)
      mf$factors <- mf$lambda.constraints <- mf$burnin <- mf$mcmc <- NULL
      mf$thin <- mf$tune <- mf$verbose <- mf$seed <- NULL
      mf$lambda.start <- mf$l0 <- mf$L0 <- NULL
      mf$store.lambda <- mf$store.scores <- mf$drop.constantvars <- NULL
      mf$... <- NULL
      mf$drop.unused.levels <- TRUE
      mf[[1]] <- as.name("model.frame")
      mf$na.action <- 'na.pass'
      mf <- eval(mf, sys.frame(sys.parent()))
      attributes(mt)$intercept <- 0
      Xterm.length <- length(attr(mt, "variables"))
      X <- subset(mf,
                  select=as.character(attr(mt, "variables"))[2:Xterm.length])
      if (drop.constantvars==TRUE){
        x.col.var <- apply(X, 2, var, na.rm=TRUE)
        X <- X[,x.col.var!=0]
      }
      N <- nrow(X)	        # number of observations      
      K <- ncol(X)              # number of manifest variables
      ncat <- matrix(NA, K, 1)  # vector of number of categ. in each man. var. 
      for (i in 1:K){
        X[,i] <- factor(X[,i], ordered=TRUE)
        ncat[i] <- nlevels(X[,i])      
        X[,i] <- as.integer(X[,i])
        X[is.na(X[,i]), i] <- -999
      }
      X <- as.matrix(X)
      xvars <- dimnames(X)[[2]] # X variable names
      xobs <- dimnames(X)[[1]]  # observation names
    }

    ## take care of the case where X has no row names
    if (is.null(xobs)){
      xobs <- 1:N
    }
    
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)
    
    ## setup constraints on Lambda
    holder <- build.factor.constraints(lambda.constraints, X, K, factors+1)
    Lambda.eq.constraints <- holder[[1]]
    Lambda.ineq.constraints <- holder[[2]]
    X.names <- holder[[3]]
  
    ## setup prior on Lambda
    holder <- form.factload.norm.prior(l0, L0, K, factors+1, X.names)
    Lambda.prior.mean <- holder[[1]]
    Lambda.prior.prec <- holder[[2]]
    if (case.switch==2){# if MCMCirtKD make it a prior on the diff. param.
      Lambda.prior.mean[,1] <- Lambda.prior.mean[,1] * -1 
    }
    
    # seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## Starting values for Lambda
    Lambda <- matrix(0, K, factors+1)
    if (is.na(lambda.start)){# sets Lambda to equality constraints & 0s
      for (i in 1:K){
        for (j in 1:(factors+1)){
          if (Lambda.eq.constraints[i,j]==-999){
            if(Lambda.ineq.constraints[i,j]==0){
              if (j==1){
                if (ncat[i] == 2){
                  probit.out <- glm(as.factor(X[X[,i]!=-999,i])~1,
                                    family=binomial(link="probit"))
                  probit.beta <- coef(probit.out)
                  Lambda[i,j] <- probit.beta[1]
                }
                if (ncat[i] > 2){
                  polr.out <- polr(ordered(X[X[,i]!=-999,i])~1)
                  Lambda[i,j] <- -polr.out$zeta[1]*.588
                }
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
    else if (is.matrix(lambda.start)){
      if (nrow(lambda.start)==K && ncol(lambda.start)==(factors+1))
        Lambda  <- lambda.start
      else {
        cat("Starting values not of correct size for model specification.\n")
        stop("Please respecify and call ", echo.name, "() again\n")
      }
    }
    else if (length(lambda.start)==1 && is.numeric(lambda.start)){
      Lambda <- matrix(lambda.start, K, factors+1)
      for (i in 1:K){
        for (j in 1:(factors+1)){
          if (Lambda.eq.constraints[i,j] != -999)
            Lambda[i,j] <- Lambda.eq.constraints[i,j]
        }
      }    
    }
    else {
      cat("Starting values neither NA, matrix, nor scalar.\n")
      stop("Please respecify and call ", echo.name, "() again\n")
    }

    ## check MH tuning parameter
    if (is.na(tune)){
      tune <- matrix(NA, K, 1)
      for (i in 1:K){
        tune[i] <- 0.05/ncat[i]
      }
    }
    else if (is.double(tune)){
      tune <- matrix(tune/ncat, K, 1)
    }
    if(min(tune) < 0) {
      cat("Tuning parameter is negative.\n")
      stop("Please respecify and call ", echo.name, "() again\n")
    }

    ## starting values for gamma (note: not changeable by user)
    gamma <- matrix(0, max(ncat)+1, K)
    for (i in 1:K){
      if (ncat[i]<=2){
        gamma[1,i] <- -300
        gamma[2,i] <- 0
        gamma[3,i] <- 300
      }
      if(ncat[i] > 2) {
        polr.out <- polr(ordered(X[X[,i]!=-999,i])~1)
        gamma[1,i] <- -300
        gamma[2,i] <- 0
        gamma[3:ncat[i],i] <- (polr.out$zeta[2:(ncat[i]-1)] -
                               polr.out$zeta[1])*.588
        gamma[ncat[i]+1,i] <- 300
      }
    }

    ## define holder for posterior sample
    if (store.scores == FALSE && store.lambda == FALSE){
      sample <- matrix(data=0, mcmc/thin, length(gamma))
    }
    else if (store.scores == TRUE && store.lambda == FALSE){
      sample <- matrix(data=0, mcmc/thin, (factors+1)*N + length(gamma))
    }
    else if(store.scores == FALSE && store.lambda == TRUE) {
      sample <- matrix(data=0, mcmc/thin, K*(factors+1)+length(gamma))
    }
    else { # store.scores==TRUE && store.lambda==TRUE
      sample <- matrix(data=0, mcmc/thin, K*(factors+1)+(factors+1)*N +
                       length(gamma))
    }


    accepts <- matrix(0, K, 1)
    
    ## Call the C++ code to do the real work
    posterior <- .C("ordfactanalpost",
                    samdata = as.double(sample),
                    samrow = as.integer(nrow(sample)),
                    samcol = as.integer(ncol(sample)),
                    X = as.integer(X),
                    Xrow = as.integer(nrow(X)),
                    Xcol = as.integer(ncol(X)),
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    tune = as.double(tune),
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                    verbose = as.integer(verbose),
                    Lambda = as.double(Lambda),
                    Lambdarow = as.integer(nrow(Lambda)),
                    Lambdacol = as.integer(ncol(Lambda)),
                    gamma = as.double(gamma),
                    gammarow = as.integer(nrow(gamma)),
                    gammacol = as.integer(ncol(gamma)),
                    ncat = as.integer(ncat),
                    ncatrow = as.integer(nrow(ncat)),
                    ncatcol = as.integer(ncol(ncat)),
                    Lameq = as.double(Lambda.eq.constraints),
                    Lameqrow = as.integer(nrow(Lambda.eq.constraints)),
                    Lameqcol = as.integer(ncol(Lambda.ineq.constraints)),
                    Lamineq = as.double(Lambda.ineq.constraints),
                    Lamineqrow = as.integer(nrow(Lambda.ineq.constraints)),
                    Lamineqcol = as.integer(ncol(Lambda.ineq.constraints)),
                    Lampmean = as.double(Lambda.prior.mean),
                    Lampmeanrow = as.integer(nrow(Lambda.prior.mean)),
                    Lampmeancol = as.integer(ncol(Lambda.prior.prec)),
                    Lampprec = as.double(Lambda.prior.prec),
                    Lampprecrow = as.integer(nrow(Lambda.prior.prec)),
                    Lamppreccol = as.integer(ncol(Lambda.prior.prec)),
                    storelambda = as.integer(store.lambda),
                    storescores = as.integer(store.scores),
                    accepts = as.integer(accepts),
                    acceptsrow = as.integer(nrow(accepts)),
                    acceptscol = as.integer(ncol(accepts)),
                    outswitch = as.integer(case.switch),
                    PACKAGE="MCMCpack"
                    )
    if(case.switch==1) {
      accepts <- matrix(posterior$accepts, posterior$acceptsrow,
                        posterior$acceptscol, byrow=FALSE)
      rownames(accepts) <- X.names
      colnames(accepts) <- ""
      cat("\n\nAcceptance rates:\n")
      print(t(accepts) / (posterior$burnin+posterior$mcmc), digits=2,
            width=6)      
    }
    
    # put together matrix and build MCMC object to return
    sample <- matrix(posterior$samdata, posterior$samrow, posterior$samcol,
                     byrow=FALSE)
    output <- mcmc(data=sample,start=1, end=mcmc, thin=thin)
    
    par.names <- NULL
    if (store.lambda==TRUE){
      if(case.switch==1) {
      Lambda.names <- paste(paste("Lambda",
                                  rep(X.names,
                                      each=(factors+1)), sep=""),
                            rep(1:(factors+1),K), sep=".")
      }
      if(case.switch==2) {
        alpha.hold <- paste("alpha", X.names, sep=".")
        beta.hold <- paste("beta", X.names, sep = ".")
        beta.hold <- rep(beta.hold, factors, each=factors)
        beta.hold <- paste(beta.hold, 1:factors, sep=".")
                
        Lambda.names <- t(cbind(matrix(alpha.hold, K, 1), 
                                matrix(beta.hold,K,factors,byrow=TRUE)))  
        dim(Lambda.names) <- NULL
      }
      par.names <- c(par.names, Lambda.names)
    }
    
    gamma.names <- paste(paste("gamma",
                               rep(0:(nrow(gamma)-1),
                                   each=K), sep=""),
                         rep(X.names,  nrow(gamma)), sep=".")
    par.names <- c(par.names, gamma.names)
    
    if (store.scores==TRUE){
      if(case.switch==1) {
      phi.names <- paste(paste("phi",
                               rep(xobs, each=(factors+1)), sep="."),
                         rep(1:(factors+1),(factors+1)), sep=".")
      par.names <- c(par.names, phi.names)
      }
      if(case.switch==2) {
      phi.names <- paste(paste("theta",
                               rep(xobs, each=(factors+1)), sep="."),
                         rep(0:factors,(factors+1)), sep=".")
      par.names <- c(par.names, phi.names)      
      
      }
    }
    
    varnames(output) <- par.names

    # get rid of columns for constrained parameters
    output.df <- as.data.frame(as.matrix(output))
    output.var <- diag(var(output.df))
    output.df <- output.df[,output.var != 0]
    output <- mcmc(as.matrix(output.df), start=burnin+1, end=burnin+mcmc,
                   thin=thin)
    
    # add constraint info so this isn't lost
    attr(output, "constraints") <- lambda.constraints
    attr(output, "n.manifest") <- K
    attr(output, "n.factors") <- factors
    attr(output, "accept.rates") <- t(accepts) / (posterior$burnin+posterior$mcmc)
    if(case.switch==1) {
      attr(output,"title") <-
        "MCMCpack Ordinal Data Factor Analysis Posterior Sample"
    }
    if(case.switch==2) {
      attr(output,"title") <-
        "MCMCpack K-Dimensional Item Response Theory Model Posterior Sample"
    }
    return(output)
    
  }

