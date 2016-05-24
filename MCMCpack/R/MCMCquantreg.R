"MCMCquantreg" <-
  function(formula, data=NULL, tau=0.5, burnin = 1000, mcmc = 10000,
           thin=1, verbose = 0, seed = sample(1:1000000,1), beta.start = NA,
           b0 = 0, B0 = 0, 
           ...) {
    
    ## checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)

    cl <- match.call()
    if (tau<=0 || tau>=1){
	stop("tau must be in (0,1). \nPlease respecify and call again.\n")
	}
    
    ## seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## form response and model matrices
    holder <- parse.formula(formula, data=data)
    Y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]    
    K <- ncol(X)  # number of covariates
        
    ## starting values and priors
    ols.fit <- lm(formula, data=data)
    defaults <- matrix(coef(ols.fit),K,1)
    defaults[1] <- defaults[1]+summary(ols.fit)$sigma*qnorm(tau)
    beta.start <- coef.start(beta.start, K, formula, family=gaussian, data, defaults=defaults)
    mvn.prior <- form.mvn.prior(b0, B0, K)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]
    

    B0.eigenvalues <- eigen(B0)$values
    if (min(B0.eigenvalues) < 0){
      stop("B0 is not positive semi-definite.\nPlease respecify and call again.\n")
    }
    

    
    ## define holder for posterior sample
    sample <- matrix(data=0, mcmc/thin, K)
    posterior <- NULL 
    
    ## call C++ code to draw samples
    auto.Scythe.call(output.object="posterior", cc.fun.name="MCMCquantreg", 
                     sample.nonconst=sample, tau=as.double(tau), Y=Y, X=X, 
                     burnin=as.integer(burnin),
                     mcmc=as.integer(mcmc), thin=as.integer(thin),
                     lecuyer=as.integer(lecuyer), 
                     seedarray=as.integer(seed.array),
                     lecuyerstream=as.integer(lecuyer.stream),
                     verbose=as.integer(verbose), betastart=beta.start,
                     b0=b0, B0=B0, package="MCMCpack")

    
## pull together matrix and build MCMC object to return
    output <- form.mcmc.object(posterior,
                               names=xnames,
                               title="MCMCquantreg Posterior Sample",
                               y=Y, call=cl)
    
    return(output)
}
