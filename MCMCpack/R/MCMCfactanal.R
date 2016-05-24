##########################################################################
## sample from the posterior distribution of a factor analysis model
## model in R using linked C++ code in Scythe.
##
## The model is:
##
## x_i = \Lambda \phi_i + \epsilon_i,   \epsilon_i \sim N(0, \Psi)
##
## where \Psi is diagonal and the priors are:
##
## \lambda_{ij} \sim N(l_{ij}, L^{-1}_{ij})
## \phi_i \sim N(0,I)
## \psi^2_{jj} \sim IG(a0_j/2, b0_j/2)
##
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
## May 7, 2003
## revised to accomodate new spec 7/5/2004 KQ
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

"MCMCfactanal" <-
  function(x, factors, lambda.constraints=list(),
           data=NULL, burnin = 1000, mcmc = 20000,
           thin=1, verbose = 0, seed = NA,
           lambda.start = NA, psi.start = NA,
           l0=0, L0=0, a0=0.001, b0=0.001,
           store.scores = FALSE, std.var=TRUE, ... ) {

    ## check for an offset
    check.offset(list(...))
    
    ## get data matrix and associated constants
    if (is.matrix(x)){
      X <- x
      xvars <- dimnames(X)[[2]]
      xobs <- dimnames(X)[[1]]
      N <- nrow(X)
      K <- ncol(X)     
    }
    else {
      holder <- parse.formula(formula=x, data=data,
                              intercept=FALSE, justX=TRUE)
      X <- holder[[2]]
      xvars <- holder[[3]]
      xobs <- holder[[4]]
      N <- nrow(X)
      K <- ncol(X)
    }

    ## standardize X
    if (std.var){
      for (i in 1:K){
        X[,i] <- (X[,i]-mean(X[,i]))/sd(X[,i])
      }
    }
    else{
      for (i in 1:K){
        X[,i] <- X[,i]-mean(X[,i])
      }
    }

    ## take care of the case where X has no row names
    if (is.null(xobs)){
      xobs <- 1:N
    }

    ## check mcmc parameters
    check.mcmc.parameters(burnin, mcmc, thin)

    # seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## setup constraints on Lambda
    holder <- build.factor.constraints(lambda.constraints, X, K, factors)
    Lambda.eq.constraints <- holder[[1]]
    Lambda.ineq.constraints <- holder[[2]]
    X.names <- holder[[3]]
  
    ## setup prior on Lambda
    holder <- form.factload.norm.prior(l0, L0, K, factors, X.names)
    Lambda.prior.mean <- holder[[1]]
    Lambda.prior.prec <- holder[[2]]

    ## setup and check prior on Psi
    holder <- form.ig.diagmat.prior(a0, b0, K)
    a0 <- holder[[1]]
    b0 <- holder[[2]]
    
    ## starting values for Lambda
    Lambda <- factload.start(lambda.start, K, factors,
                             Lambda.eq.constraints, Lambda.ineq.constraints)
    
    ## starting values for Psi
    Psi <- factuniqueness.start(psi.start, X)
    

    ## define holder for posterior sample
    if(store.scores == FALSE) {
      sample <- matrix(data=0, mcmc/thin, K*factors+K)
    }
    else {
      sample <- matrix(data=0, mcmc/thin, K*factors+K+factors*N)
    }
    posterior <- NULL 
    
    ## call C++ code to do the sampling
    auto.Scythe.call(output.object="posterior", cc.fun.name="MCMCfactanal",
                     sample.nonconst=sample, X=X, burnin=as.integer(burnin),
                     mcmc=as.integer(mcmc), thin=as.integer(thin), 
                     lecuyer=as.integer(lecuyer), 
                     seedarray=as.integer(seed.array),
                     lecuyerstream=as.integer(lecuyer.stream),
                     verbose=as.integer(verbose),
                     Lambda=Lambda, Psi=Psi, Lameq=Lambda.eq.constraints,
                     Lamineq=Lambda.ineq.constraints,
                     Lampmean=Lambda.prior.mean, Lampprec=Lambda.prior.prec,
                     a0=a0, b0=b0, storescores=as.integer(store.scores))

    Lambda.names <- paste(paste("Lambda",
                                rep(X.names,
                                    each=factors), sep=""),
                          rep(1:factors,K), sep="_")
    Psi.names <- paste("Psi", X.names, sep="")
    par.names <- c(Lambda.names, Psi.names)
    if (store.scores==TRUE){
      phi.names <- paste(paste("phi",
                               rep(xobs, each=factors), sep="_"),
                         rep(1:factors,factors), sep="_")
      par.names <- c(par.names, phi.names)
    }
    title <- "MCMCpack Factor Analysis Posterior Sample"
    output <- form.mcmc.object(posterior, par.names, title)
    ## get rid of columns for constrained parameters
    output.df <- as.data.frame(as.matrix(output))
    output.sd <- apply(output.df, 2, sd)
    output.df <- output.df[,output.sd != 0]
    
    output <- mcmc(as.matrix(output.df), start=burnin+1, end=mcmc+burnin,
                   thin=thin)

    ## add constraint info so this isn't lost
    attr(output, "constraints") <- lambda.constraints
    attr(output, "n.manifest") <- K
    attr(output, "n.factors") <- factors
    return(output)
  }
