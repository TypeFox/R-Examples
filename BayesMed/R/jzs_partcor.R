jzs_partcor <-
  function(V1,V2,control,
           alternative=c("two.sided","less","greater"),
           n.iter=10000,n.burnin=500,standardize=TRUE){
    
    runif(1) # defines .Random.seed
    
    # standardize variables
    if(standardize==TRUE){
      M <- (V1-mean(V1))/sd(V1)
      Y <- (V2-mean(V2))/sd(V2)
      X <- (control-mean(control))/sd(control)
    } else {
      M <- V1
      Y <- V2
      X <- control      
    }
    
    r0 <- sqrt(summary(lm(M~X))$r.squared)
    r1 <- sqrt(summary(lm(M~X+Y))$r.squared)
    p0 <- 1
    p1 <- 2
    n  <- length(X)
    
    # main function to analytically calculate the BF for partial correlation
    # see Wetzels, R. & Wagenmakers, E.-J. (2012). A default Bayesian hypothesis test for correlations and partial correlations. Psychonomic Bulletin & Review.
    jzs_partcorbf <- function(r0,r1,p0,p1,n){
      int <- function(r,n,p,g){
        a <- .5*((n-1-p)*log(1+g)-(n-1)*log(1+g*(1-r^2)))
        exp(a)*dinvgamma(g,shape=.5,scale=n/2)
      }
      bf10 <- integrate(int,lower=0,upper=Inf,r=r1,p=p1,n=n)$value/
        integrate(int,lower=0,upper=Inf,r=r0,p=p0,n=n)$value
      return(bf10)
    }
    
    BF <- jzs_partcorbf(r0,r1,p0,p1,n)
    
    ### the next part is needed to impose order restrictions
    ### for the order restrictions we need to estimate the posterior samples
    
    #==========================================================
    # load JAGS models
    #==========================================================
    
    
    jagsmodelpartialcorrelation <- 
      
      "####### Cauchy-prior on beta and tau' #######
    model
    
{
    
    for (i in 1:n)
    
{
    mu[i] <- intercept + theta[1]*x[i,1] + theta[2]*x[i,2]
    y[i]   ~ dnorm(mu[i],phi)
    
}
    
    # uninformative prior on intercept alpha, 
    # Jeffreys' prior on precision phi
    intercept ~ dnorm(0,.0001)
    phi   ~ dgamma(.0001,.0001)
    #phi   ~ dgamma(0.0000001,0.0000001) #JAGS accepts even this
    #phi   ~ dgamma(0.01,0.01)           #WinBUGS wants this
    
    # inverse-gamma prior on g:
    g       <- 1/invg 
    a.gamma <- 1/2
    b.gamma <- n/2    
    invg     ~ dgamma(a.gamma,b.gamma)
    
    # Ntzoufras, I. (2009). Bayesian Modeling Using WinBUGS.
    # New Jersey: John Wiley & Sons, Inc. p. 167
    # calculation of the inverse matrix of V
    inverse.V <- inverse(V)
    # calculation of the elements of prior precision matrix
    for(i in 1:2)
{ 
    for (j in 1:2)
{
    prior.T[i,j] <- inverse.V[i,j] * phi/g
}
}
    # multivariate prior for the beta vector
    theta[1:2] ~ dmnorm( mu.theta, prior.T )
    for(i in 1:2) { mu.theta[i] <- 0 }
    
}
    
    # Explanation-----------------------------------------------------------------
    # Prior on g:
    # We know that g ~ inverse_gamma(1/2, n/2), with 1/2 the shape parameter and 
    # n/2 the scale parameter.
    # It follows that 1/g ~ gamma(1/2, 2/n).
    # However, BUGS/JAGS uses the *rate parameterization* 1/theta instead of the
    # scale parametrization theta. Hence we obtain, in de BUGS/JAGS rate notation:
    # 1/g ~ dgamma(1/2, n/2)
    # Also note: JAGS does not want [,] structure
    #-----------------------------------------------------------------------------
    "
    
    jags.model.file2 <- tempfile(fileext=".txt")
    write(jagsmodelpartialcorrelation,jags.model.file2)
    
    #==========================================================
    # BF FOR PARTIAL CORRELATION (MY|X)
    #==========================================================
    
    x <- cbind(X,M)
    y <- Y
    
    V <- solve(t(x)%*%x) #NB I switched to the notation from Ntzoufras, p. 167
    
    jags.data   <- list("n", "x", "y", "V")
    jags.params <- c("theta")
    jags.inits  <-  list(
      list(theta = c(0.0,0.3)),  #chain 1 starting value
      list(theta = c(0.3, 0.0)), #chain 2 starting value
      list(theta = c(-.15,.15))) #chain 3 starting value
    
    jagssamples <- jags(data=jags.data, inits=jags.inits, jags.params, 
                        n.chains=3, n.iter=n.iter, DIC=T,
                        n.burnin=n.burnin, n.thin=1, model.file=jags.model.file2)
    
    beta <- jagssamples$BUGSoutput$sims.list$theta[,2]
    
    #-------------------------------------------------------
    
    # one-sided test?
    
    # save BF for one-tailed test
    # BF21 = 2*{proportion posterior samples of beta < 0}
    propposterior_less <- sum(beta<0)/length(beta)
    propposterior_greater <- sum(beta>0)/length(beta)
    
    # posterior proportion cannot be zero, because this renders a BF of zero
    # none of the samples of the parameter follow the restriction
    # ergo: the posterior proportion is smaller than 1/length(parameter)
    
    if(propposterior_less==0){
      propposterior_less <- 1/length(beta)
    }
    
    
    if(propposterior_greater==0){
      propposterior_greater <- 1/length(beta)
    }
    
    BF21_less <- 2*propposterior_less
    BF21_greater <- 2*propposterior_greater
    
    if(alternative[1]=="less"){
      # BF10 = p(D|b~cauchy(0,1))/p(D|b=0)
      BF10 <- BF
      
      # BF21 = p(D|b~cauchy-(0,1))/p(D|b~cauchy(0,1))
      # BF21 = 2*{proportion posterior samples of beta < 0}
      BF21 <- BF21_less
      
      BF <- BF10*BF21
      
    } else if(alternative[1]=="greater"){
      # BF10 = p(D|b~cauchy(0,1))/p(D|b=0)
      BF10 <- BF
      
      # BF21 = p(D|b~cauchy+(0,1))/p(D|b~cauchy(0,1))
      # BF21 = 2*{proportion posterior samples of beta > 0}
      BF21 <- BF21_greater
      
      BF <- BF10*BF21
      
    }
    
    #---------------------------------------------------
    
    # convert BFs to posterior probability
    # prob cannot be exactly 1 or 0
    prob_b <- BF/(BF+1)
    
    if(prob_b == 1){
      prob_b <- prob_b - .Machine$double.eps
    }
    if(prob_b == 0){
      prob_b <- prob_b + .Machine$double.eps
    }
    
    
    #====================================================
    
    res <- list(PartCoef=mean(beta),
                BayesFactor=BF,
                PosteriorProbability=prob_b,
                beta_samples=beta,
                jagssamples=jagssamples)
    
    class(res) <- c("jzs_med","list")
    class(res$jagssamples) <- "rjags"
    class(res$beta_samples) <- "CI"
    
    return(res)
    
  }
