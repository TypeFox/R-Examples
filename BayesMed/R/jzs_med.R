jzs_med <-
  function(independent,dependent,mediator,
           alternativeA=c("two.sided","less","greater"),
           alternativeB=c("two.sided","less","greater"),
           alternativeT=c("two.sided","less","greater"),
           n.iter=10000,n.burnin=500,
           standardize=TRUE){
    
    runif(1) # defines .Random.seed
    
    # independent = vector with values for independent variable
    # dependent = vector with values for dependent variable
    # mediator = vector with values for mediating variable
    
    # sample size
    n <- length(independent) 
    
    X <- independent
    Y <- dependent
    M <- mediator
    
    if(standardize==TRUE){
      X <- (X-mean(X))/sd(X)
      Y <- (Y-mean(Y))/sd(Y)
      M <- (M-mean(M))/sd(M)
    }
    
    #==========================================================
    # RESULTS FOR PATH ALPHA 
    #==========================================================
    
    #     alternativeA <- alternativeA
    #     n.iter <- n.iter
    #     n.burnin <- n.burnin
    
    res_alpha    <- jzs_cor(X,M,alternative=alternativeA,n.iter=n.iter,
                            n.burnin=n.burnin,standardize=standardize)
    
    BFa <- res_alpha$BayesFactor
    prob_a <- res_alpha$PosteriorProbability
    alpha <- res_alpha$alpha
    jagssamplesA <- res_alpha$jagssamples
    
    #==========================================================
    
    # we chose not to use jzs_partcor() for the results of path beta and tau
    # because this would mean we'd have to fit the JAGS model three times in total
    # now we can get the information for both beta and tau out of one and the same model
    
    # function to analytically calculate the BF for partial correlation
    # see Wetzels, R. & Wagenmakers, E.-J. (2012). A default Bayesian hypothesis test for correlations and partial correlations. Psychonomic Bulletin & Review.
    jzs_partcor_basic <- function(V1,V2,control,standardize=TRUE){
      
      # standardize variables
      if(standardize==TRUE){
        V1 <- (V1-mean(V1))/sd(V1)
        V2 <- (V2-mean(V2))/sd(V2)
        control <- (control-mean(control))/sd(control)
      }
      
      r0 <- sqrt(summary(lm(V1~control))$r.squared)
      r1 <- sqrt(summary(lm(V1~control+V2))$r.squared)
      p0 <- 1
      p1 <- 2
      n  <- length(V1)
      
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
      return(BF)
    }
    
    #==========================================================
    # JAGS MODEL FOR PARTIAL CORRELATION
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
    # SAVE SAMPLES FOR PATH BETA AND TAU_ACCENT
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
    
    jagssamplesTB <- jags(data=jags.data, inits=jags.inits, jags.params, 
                          n.chains=3, n.iter=n.iter, DIC=T,
                          n.burnin=n.burnin, n.thin=1, model.file=jags.model.file2)
    
    tau_accent <- jagssamplesTB$BUGSoutput$sims.list$theta[,1]
    beta <- jagssamplesTB$BUGSoutput$sims.list$theta[,2]
    
    #==========================================================
    # RESULTS FOR PATH BETA 
    #==========================================================
    
    BFb <- jzs_partcor_basic(M,Y,control=X,standardize=standardize)
    
    # one-sided test beta?
    
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
    
    if(alternativeB[1]=="less"){
      # BF10 = p(D|b~cauchy(0,1))/p(D|b=0)
      BF10 <- BFb
      
      # BF21 = p(D|b~cauchy-(0,1))/p(D|b~cauchy(0,1))
      # BF21 = 2*{proportion posterior samples of beta < 0}
      BF21 <- BF21_less
      
      BFb <- BF10*BF21
      
    } else if(alternativeB[1]=="greater"){
      # BF10 = p(D|b~cauchy(0,1))/p(D|b=0)
      BF10 <- BFb
      
      # BF21 = p(D|b~cauchy+(0,1))/p(D|b~cauchy(0,1))
      # BF21 = 2*{proportion posterior samples of beta > 0}
      BF21 <- BF21_greater
      
      BFb <- BF10*BF21
      
    }
    
    #---------------------------------------------------
    
    # convert BFb to posterior probability
    # prob cannot be exactly 1 or 0
    prob_b <- BFb/(BFb+1)
    
    if(prob_b == 1){
      prob_b <- prob_b - .Machine$double.eps
    }
    if(prob_b == 0){
      prob_b <- prob_b + .Machine$double.eps
    }
    
    #=========================================
    # calculate evidence for mediation (EM)
    #=========================================
    
    EM <- prob_a*prob_b
    BF.EM <- EM/(1-EM)
    
    #==========================================================
    # RESULTS FOR PATH TAU_ACCENT 
    #==========================================================
    
    BFt_accent <- jzs_partcor_basic(X,Y,control=M,standardize=standardize)
    
    # one-sided test tau_accent?
    
    # save BF for one-tailed test
    # BF21 = 2*{proportion posterior samples of tau_accent < 0}
    propposterior_less <- sum(tau_accent<0)/length(tau_accent)
    propposterior_greater <- sum(tau_accent>0)/length(tau_accent)
    
    # posterior proportion cannot be zero, because this renders a BF of zero
    # none of the samples of the parameter follow the restriction
    # ergo: the posterior proportion is smaller than 1/length(parameter)
    
    if(propposterior_less==0){
      propposterior_less <- 1/length(tau_accent)
    }
    
    
    if(propposterior_greater==0){
      propposterior_greater <- 1/length(tau_accent)
    }
    
    BF21_less <- 2*propposterior_less
    BF21_greater <- 2*propposterior_greater
    
    if(alternativeT[1]=="less"){
      # BF10 = p(D|t'~cauchy(0,1))/p(D|t'=0)
      BF10 <- BFt_accent
      
      # BF21 = p(D|t'~cauchy-(0,1))/p(D|t'~cauchy(0,1))
      # BF21 = 2*{proportion posterior samples of tau_accent < 0}
      BF21 <- BF21_less
      
      BFt_accent <- BF10*BF21
      
    } else if(alternativeT[1]=="greater"){
      # BF10 = p(D|t'~cauchy(0,1))/p(D|t'=0)
      BF10 <- BFt_accent
      
      # BF21 = p(D|t'~cauchy+(0,1))/p(D|t'~cauchy(0,1))
      # BF21 = 2*{proportion posterior samples of tau_accent > 0}
      BF21 <- BF21_greater
      
      BFt_accent <- BF10*BF21
      
    }
    
    #--------------------------------------------------------
    
    # convert BFs to posterior probability
    # prob cannot be exactly 1 or 0
    prob_t_accent <- BFt_accent/(BFt_accent+1)
    
    if(prob_t_accent == 1){
      prob_t_accent <- prob_t_accent - .Machine$double.eps
    }
    if(prob_t_accent == 0){
      prob_t_accent <- prob_t_accent + .Machine$double.eps
    }
    
    #===============================================================
    
    # calculate 95% credible interval for ab
    ab <- alpha*beta
    CI <- quantile(ab,c(.025,.975))
    
    #===============================================================
    
    res <- data.frame(Estimate = c(mean(alpha),mean(beta),mean(tau_accent),mean(ab)),
                      BF = c(BFa,BFb,BFt_accent,BF.EM),
                      PostProb = c(prob_a,prob_b,prob_t_accent,EM))
    
    rownames(res) <- c("alpha","beta","tau_prime","Mediation (alpha*beta)")
    
    result <- list(main_result=res,
                   CI_ab=CI,
                   alpha_samples=alpha,
                   beta_samples=beta,
                   tau_prime_samples=tau_accent,
                   ab_samples=ab,
                   jagssamplesA=jagssamplesA,
                   jagssamplesTB=jagssamplesTB)
    
    
    class(result) <- c("jzs_med","list")
    class(result$main_result) <- c("JZSMed","data.frame")
    class(result$jagssamplesA) <- "rjags"
    class(result$jagssamplesTB) <- "rjags"
    class(result$ab_samples) <- "CI"
    class(result$alpha_samples) <- "CI"
    class(result$beta_samples) <- "CI"
    class(result$tau_prime_samples) <- "CI"
    
    return(result)
    
  }
