jzs_medSD <-
  function(independent,dependent,mediator,
           SDmethod=c("fit.st","dnorm","splinefun","logspline"),
           alternativeA=c("two.sided","less","greater"),
           alternativeB=c("two.sided","less","greater"),
           alternativeT=c("two.sided","less","greater"),
           n.iter=10000,n.burnin=500,
           standardize=TRUE){
    
    runif(1) # defines .Random.seed
    
    # independent = vector with values for independent variable
    # dependent = vector with values for dependent variable
    # mediator = vector with values for mediating variable
    
    X <- independent
    Y <- dependent
    M <- mediator
    
    if(standardize==TRUE){
      X <- (X-mean(X))/sd(X)
      Y <- (Y-mean(Y))/sd(Y)
      M <- (M-mean(M))/sd(M)
    }
      
    # sample size
    n <- length(independent)
    
    #==========================================================
    # load JAGS models
    #==========================================================
    
    jagsmodelcorrelation <- 
      
      "####### Cauchy-prior on single beta #######
model
    
{
    
    for (i in 1:n)
    
{
    mu[i] <- intercept + alpha*x[i]
    y[i]   ~ dnorm(mu[i],phi)
    
}
    
    # uninformative prior on the intercept intercept, 
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
    
    
    # g-prior on beta:
    vari <- (g/phi) * invSigma 
    prec <- 1/vari
    alpha    ~ dnorm(0, prec)
}
    
    # Explanation------------------------------------------------------------------ 
    # Prior on g:
    # We know that g ~ inverse_gamma(1/2, n/2), with 1/2 the shape
    # parameter and n/2 the scale parameter.
    # It follows that 1/g ~ gamma(1/2, 2/n).
    # However, BUGS/JAGS uses the *rate parameterization* 1/theta instead of the
    # scale parametrization theta. Hence we obtain, in de BUGS/JAGS rate notation:
    # 1/g ~ dgamma(1/2, n/2)
    #------------------------------------------------------------------------------
    "
    jags.model.file1 <- tempfile(fileext=".txt")
    write(jagsmodelcorrelation,jags.model.file1)
    
    #============================================================================
    
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
    # BF FOR PATH alpha: CORRELATION X-M
    #==========================================================
    
    x <- X
    y <- M
    
    invSigma <- solve(t(x)%*%x)
    
    jags.data   <- list("n", "x", "y", "invSigma")
    jags.params <- c("alpha", "g")
    jags.inits  <-  list(
      list(alpha = 0.0), #chain 1 starting value
      list(alpha = -0.3), #chain 2 starting value
      list(alpha = 0.3)) #chain 3 starting value
    
    jagssamplesA <- jags(data=jags.data, inits=jags.inits, jags.params, 
                       n.chains=3, n.iter=n.iter, DIC=T,
                       n.burnin=n.burnin, n.thin=1, model.file=jags.model.file1)
    
    # estimate the posterior regression coefficient and scaling factor g
    alpha <- jagssamplesA$BUGSoutput$sims.list$alpha[,1]
    g  <- jagssamplesA$BUGSoutput$sims.list$g
    
    #------------------------------------------------------------------
    
    if(SDmethod[1]=="fit.st"){
      
      mydt <- function(x, m, s, df) dt((x-m)/s, df)/s
      
      foo <- try({
        fit.t1 <- fit.st(alpha)
        nuA    <- as.numeric(fit.t1$par.ests[1]) #degrees of freedom
        muA    <- as.numeric(fit.t1$par.ests[2]) 
        sigmaA <- abs(as.numeric(fit.t1$par.ests[3])) # This is a hack -- with high n occasionally
        # sigma switches sign. 
      }) 
      
      if(!("try-error"%in%class(foo))){
        
        # BAYES FACTOR ALPHA
        BFa <- 1/(mydt(0,muA,sigmaA,nuA)/dcauchy(0))
        
        # save BF for one-tailed test
        # BF21 = 2*{proportion posterior samples of alpha < 0}
        BF21a_less <- 2*pt((0-muA)/sigmaA,nuA,lower.tail=TRUE)/sigmaA
        BF21a_greater <- 2*pt((0-muA)/sigmaA,nuA,lower.tail=FALSE)/sigmaA
        
      } else {
        
        warning("fit.st did not converge, alternative optimization method was used.","\n")
        
        mydt2 <- function(pars){
          
          mA <- pars[1]
          sA <- abs(pars[2])  # no negative standard deviation
          dfA <- abs(pars[3]) # no negative degrees of freedom
          
          -2*sum(dt((alpha-mA)/sA, dfA,log=TRUE)-log(sA))
        }
        
        res <- optim(c(mean(alpha),sd(alpha),20),mydt2)$par
        
        mA <- res[1]
        sA <- res[2]
        dfA <- res[3]
        
        
        # ALTERNATIVE BAYES FACTOR ALPHA
        BFa <- 1/(mydt2(0,mA,sA,dfA)/dcauchy(0))
        
      }
      
      #-------------------------
      
    } else if(SDmethod[1]=="dnorm"){
      BFa <- 1/(dnorm(0,mean(alpha),sd(alpha))/dcauchy(0))
       
      #-------------------------
      
    } else if(SDmethod[1]=="splinefun"){
      f <- splinefun(density(alpha))
      BFa <- 1/(f(0)/dcauchy(0))
            
      #-------------------------
      
    } else if (SDmethod[1]=="logspline"){
      fit.posterior <- logspline(alpha)
      posterior.pp  <- dlogspline(0, fit.posterior) # this gives the pdf at point b2 = 0
      prior.pp      <- dcauchy(0)                   # height of prior at b2 = 0
      BFa           <- prior.pp/posterior.pp
      
    } 
    
    #--------------------------------------------------------
    
    # one-sided test?
    
    # save BF for one-tailed test
    # BF21 = 2*{proportion posterior samples of alpha < 0}
    propposterior_less <- sum(alpha<0)/length(alpha)
    propposterior_greater <- sum(alpha>0)/length(alpha)
    
    # posterior proportion cannot be zero, because this renders a BF of zero
    # none of the samples of the parameter follow the restriction
    # ergo: the posterior proportion is smaller than 1/length(parameter)
    
    if(propposterior_less==0){
      propposterior_less <- 1/length(alpha)
    }
    
    if(propposterior_greater==0){
      propposterior_greater <- 1/length(alpha)
    }
    
    BF21a_less <- 2*propposterior_less
    BF21a_greater <- 2*propposterior_greater
    
    if(alternativeA[1]=="less"){
      # BF10 = p(D|a~cauchy(0,1))/p(D|a=0)
      BF10 <- BFa
      
      # BF21 = p(D|a~cauchy-(0,1))/p(D|a~cauchy(0,1))
      # BF21 = 2*{proportion posterior samples of alpha < 0}
      BF21 <- BF21a_less
      
      BFa <- BF10*BF21
      
    } else if(alternativeA[1]=="greater"){
      # BF10 = p(D|a~cauchy(0,1))/p(D|a=0)
      BF10 <- BFa
      
      # BF21 = p(D|a~cauchy+(0,1))/p(D|a~cauchy(0,1))
      # BF21 = 2*{proportion posterior samples of alpha > 0}
      BF21 <- BF21a_greater
      
      BFa <- BF10*BF21
      
    }
    
    #--------------------------------------------------------
    
    # convert BFs to posterior probability
    # prob cannot be exactly 1 or 0
    prob_a <- BFa/(BFa+1)
    
    if(prob_a == 1){
      prob_a <- prob_a - .Machine$double.eps
    }
    if(prob_a == 0){
      prob_a <- prob_a + .Machine$double.eps
    }
    
    #==========================================================
    # BF FOR PATH beta: PARTIAL CORRELATION MY|X
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
    
    beta <- jagssamplesTB$BUGSoutput$sims.list$theta[,2]
    
    #------------------------------------------------------------------
    
    if(SDmethod[1]=="fit.st"){
      
      bar <- try({
        fit.t2 <- fit.st(beta)
        nuB    <- as.numeric(fit.t2$par.ests[1]) #degrees of freedom
        muB    <- as.numeric(fit.t2$par.ests[2]) 
        sigmaB <- abs(as.numeric(fit.t2$par.ests[3])) # This is a hack -- with high n occasionally
        # sigma switches sign. 
      })
      
      if(!("try-error"%in%class(bar))){
        
        # BAYES FACTOR BETA
        BFb <- 1/(mydt(0,muB,sigmaB,nuB)/dcauchy(0))
        
      } else {
        
        warning("fit.st did not converge, alternative optimization method was used.","\n")
        
        mydt2 <- function(pars){
          
          mB <- pars[1]
          sB <- abs(pars[2])  # no negative standard deviation
          dfB <- abs(pars[3]) # no negative degrees of freedom
          
          -2*sum(dt((beta-mB)/sB, dfB,log=TRUE)-log(sB))
        }
        
        res <- optim(c(mean(beta),sd(beta),20),mydt2)$par
        
        mB <- res[1]
        sB <- res[2]
        dfB <- res[3]
        
        # ALTERNATIVE BAYES FACTOR BETA
        BFb <- 1/(mydt2(0,mB,sB,dfB)/dcauchy(0))
      }
      
      #-------------------------
      
    } else if(SDmethod[1]=="dnorm"){
      
      BFb <- 1/(dnorm(0,mean(beta),sd(beta))/dcauchy(0)) 
      
      #-------------------------
      
    } else if(SDmethod[1]=="splinefun"){
      f <- splinefun(density(beta))
      BFb <- 1/(f(0)/dcauchy(0))
       
      #-------------------------
      
    } else if (SDmethod[1]=="logspline"){
      fit.posterior <- logspline(beta)
      posterior.pp  <- dlogspline(0, fit.posterior) # this gives the pdf at point b2 = 0
      prior.pp      <- dcauchy(0)                   # height of prior at b2 = 0
      BFb           <- prior.pp/posterior.pp
      
    }
    
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
    
    BF21b_less <- 2*propposterior_less
    BF21b_greater <- 2*propposterior_greater
    
    if(alternativeB[1]=="less"){
      # BF10 = p(D|b~cauchy(0,1))/p(D|b=0)
      BF10 <- BFb
      
      # BF21 = p(D|b~cauchy-(0,1))/p(D|b~cauchy(0,1))
      # BF21 = 2*{proportion posterior samples of beta < 0}
      BF21 <- BF21b_less
      
      BFb <- BF10*BF21
      
    } else if(alternativeB[1]=="greater"){
      # BF10 = p(D|b~cauchy(0,1))/p(D|b=0)
      BF10 <- BFb
      
      # BF21 = p(D|b~cauchy+(0,1))/p(D|b~cauchy(0,1))
      # BF21 = 2*{proportion posterior samples of beta > 0}
      BF21 <- BF21b_greater
      
      BFb <- BF10*BF21
      
    }
    
    #---------------------------------------------------
    
    # convert BFs to posterior probability
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
    
    
    #=========================================
    # FULL OR PARTIAL MEDIATION
    #=========================================
    
    tau_accent <- jagssamplesTB$BUGSoutput$sims.list$theta[,1]
    
    #---------------------------------------------------
    
    if(SDmethod[1]=="fit.st"){
      
      baz <- try({
        fit.t3 <- fit.st(tau_accent)
        nuT    <- as.numeric(fit.t3$par.ests[1]) #degrees of freedom
        muT    <- as.numeric(fit.t3$par.ests[2]) 
        sigmaT <- abs(as.numeric(fit.t3$par.ests[3])) # This is a hack -- with high n occasionally
        # sigma switches sign. 
      })
      
      if(!("try-error"%in%class(baz))){
        
        # BAYES FACTOR TAU
        BFt_accent <- 1/(mydt(0,muT,sigmaT,nuT)/dcauchy(0))
              
      } else {
        
        warning("fit.st did not converge. Alternative optimization method was used.","\n")
        
        mydt2 <- function(pars){
          
          mT <- pars[1]
          sT <- abs(pars[2])  # no negative standard deviation
          dfT <- abs(pars[3]) # no negative degrees of freedom
          
          -2*sum(dt((tau_accent-mT)/sT, dfT,log=TRUE)-log(sT))
        }
        
        res <- optim(c(mean(tau_accent),sd(tau_accent),20),mydt2)$par
        
        mT <- res[1]
        sT <- res[2]
        dfT <- res[3]
        
        # ALTERNATIVE BAYES FACTOR TAU
        BFt_accent <- 1/(mydt2(0,mT,sT,dfT)/dcauchy(0))
       
      }
      
      #-------------------------
      
    } else if(SDmethod[1]=="dnorm"){
      BFt_accent <- 1/(dnorm(0,mean(tau_accent),sd(tau_accent))/dcauchy(0))  
      
      #-------------------------
      
    } else if(SDmethod[1]=="splinefun"){
      f <- splinefun(density(tau_accent))
      BFt_accent <- 1/(f(0)/dcauchy(0))
     
      #-------------------------
      
    } else if (SDmethod[1]=="logspline"){
      fit.posterior <- logspline(tau_accent)
      posterior.pp  <- dlogspline(0, fit.posterior) # this gives the pdf at point b2 = 0
      prior.pp      <- dcauchy(0)                   # height of prior at b2 = 0
      BFt_accent <- prior.pp/posterior.pp
      
    }
    
    #------------------------------------------------------------
    
    # one-sided test?
    
    # save BF for one-tailed test
    # BF21 = 2*{proportion posterior samples of tau < 0}
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
    
    BF21t_less <- 2*propposterior_less
    BF21t_greater <- 2*propposterior_greater
    
    
    if(alternativeT[1]=="less"){
      # BF10 = p(D|t~cauchy(0,1))/p(D|t=0)
      BF10 <- BFt_accent
      
      # BF21 = p(D|t~cauchy-(0,1))/p(D|t~cauchy(0,1))
      # BF21 = 2*{proportion posterior samples of tau_accent < 0}
      BF21 <- BF21t_less
      
      BFt_accent <- BF10*BF21
      
    } else if(alternativeT[1]=="greater"){
      # BF10 = p(D|t~cauchy(0,1))/p(D|t=0)
      BF10 <- BFt_accent
      
      # BF21 = p(D|t~cauchy+(0,1))/p(D|t~cauchy(0,1))
      # BF21 = 2*{proportion posterior samples of tau_accent > 0}
      BF21 <- BF21t_greater
      
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
    
    if(BFa<0){
      BFa <- NA
      warning("Negative Bayes factor: try other SDmethod","\n")
    }
    
    if(BFb<0){
      BFb <- NA
      warning("Negative Bayes factor: try other SDmethod","\n")
    }
    
    if(BFt_accent<0){
      warning("Negative Bayes factor: try other SDmethod","\n")
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
