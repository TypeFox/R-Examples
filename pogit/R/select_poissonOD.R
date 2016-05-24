#### Bayesian variable selection via Dirac spike and Student-t/ or normal slab 
#### for Poisson regression models with observation-specific random intercept 
#### to accommodate overdispersion in count data

#### (invisible function)
#### last change: 2016/03/25

####  ------------------------------------------------------------------------ #
select_poissonOD <- function(y, X, offset, H = NULL, mcomp, compmix.pois = NULL, 
                             cm1, model, prior, mcmc, param, imc){
  
  # linear predicor in Poisson model
  muP <- X%*%param$beta
  if (model$ri == 1){
    ranEff <- param$btilde[model$Zp]
    linp <- muP + ranEff*param$theta 
  } else linp <- muP  
  
  n <- length(y)
  
  #### Step A --- data augmentation for the Poisson model 
  lambda  <- exp(linp)
  
  ## A.1 sample the inter-arrival times of assumed Poisson process in [0,1]
  
  # IAMS - mixture components
  if (model$family == "pogit"){
    compmix.pois <- get_mixcomp(y, mcomp)
    #compmix.pois <- do.call(get_mixcomp, list(y, mcomp))  
  } 
  
  augPois <- iams1_poisson(y, offset*lambda, compmix.pois)
  tau1    <- augPois$t1
  tau2    <- augPois$t2
  
  logMu   <- linp + log(offset)
  logMugz <- logMu[compmix.pois$igz]
  
  ## A.2 --- sample the component indicators
  R <- iams2_poisson(n, tau1, tau2, logMu, logMugz, cm1, compmix.pois)
  
  # mixture component means and variances
  m1 <- cm1$comp$m[R[1:n]]
  m2 <- compmix.pois$my[cbind(seq_len(compmix.pois$ngz), R[-(1:n)])]
  mR <- as.matrix(c(m1,m2), (n + compmix.pois$ngz))
  v1 <- cm1$comp$v[R[1:n]]
  v2 <- compmix.pois$vy[cbind(seq_len(compmix.pois$ngz), R[-(1:n)])]
  invSig <- 1/sqrt(c(v1,v2))
  
  # stacking and standardizing
  tauS <- c(tau1, tau2)
  offsetS <- c(offset, offset[compmix.pois$igz])
  xS <- rbind(X, X[compmix.pois$igz, , drop = FALSE])
  yS <- (-log(tauS) - mR - log(offsetS))*invSig
  
  if (model$ri==0){
    Xall <- xS*kronecker(matrix(1, 1, model$d + 1), invSig)
  } else {
    Xall <- cbind(xS, c(ranEff, ranEff[compmix.pois$igz]))*kronecker(matrix(1, 1, model$d + model$ri + 1), invSig)
  }
  
  # inverse prior variance of regression effects (updated)
  invA0 <- diag(c(prior$invM0, 1/param$psi), nrow = model$d + model$ri + 1)
  
  
  #### Step B --- starts Bayesian variable selection
  if (imc > mcmc$startsel && sum(sum(!model$deltafix) + sum(!model$gammafix)) > 0){
    
    ## (1) update mixture weights
    incfix <- sum(param$delta==1)  
    omega  <- rbeta(1, prior$w['wa0'] + incfix, prior$w['wb0'] + model$d - incfix)
    
    if (model$ri==1){
      incran <- sum(param$gamma==1) 
      pi     <- rbeta(1, prior$pi['pa0'] + incran, prior$pi['pb0'] + model$ri - incran)
    } else {
      pi <- NULL
    }     
    
    ## (2) sample the indicators, the regression coefficients and the scale parameters
    ## --- (i) sample the indicators delta_{beta,j}, gamma_{beta} for the slab component
    
    ##==TODO generalize draw_indices for logit, pogit, poisson!!!
    indic <- draw_indicators(yS, Xall, param$delta, param$gamma, omega, pi, model, prior, invA0)
    delta  <- indic$deltanew 
    pdelta <- indic$pdeltanew
    gamma  <- indic$gammanew
    pgamma <- indic$pgammanew
  } else {
    delta  <- param$delta
    pdelta <- param$pdelta
    gamma  <- param$gamma
    pgamma <- param$pgamma
    omega  <- param$omega
    pi     <- param$pi
  }
  
  ## --- (ii_A) sample the (selected) regression effects 
  if (model$ri == 0){
    index <- c(1, which(delta == 1) + 1) 
  } else {
    index <- c(1, which(c(delta, gamma) == 1) + 1)  
  }
  
  Zsel      <- Xall[, index, drop=FALSE]  # Z*=[1, W^delta,atilde]*sqrt(Sigma^-1)
  dsel      <- length(index)
  invA0_sel <- invA0[index, index, drop=FALSE]
  a0_sel    <- prior$a0[index,,drop=FALSE]
  
  AP    <- solve(invA0_sel + t(Zsel)%*%Zsel)  # A = (A0^-1 + (Z*)'Sigma^-1 Z*)
  aP    <- AP%*%(invA0_sel%*%a0_sel + t(Zsel)%*%yS) # a = A(A0^-1*a0 + (Z*)'Sigma^-1*y
  zetaP <- t(chol(AP))%*%matrix(rnorm(dsel), dsel, 1) + aP 
  
  val1        <- matrix(0, 1, model$d + model$ri + 1) 
  val1[index] <- t(zetaP)
  mu_beta  <- val1[1]
  
  if (model$d > 0){
    beta  <- val1[2:(model$d+1)]
  } else {
    beta  <- matrix(0, 1, model$d)   
  }
  
  muP <- X%*%c(mu_beta, beta)
  
  ## --- (ii_B) sample the random intercepts
  if (model$ri == 1){
    theta <- val1[((model$d + 1) + 1):(model$d + model$ri + 1)]   
    
    yh <- (-log(tauS) - mR - log(offsetS) - c(muP, muP[compmix.pois$igz]))
    v2inv <- yh2 <- rep(0, n)
    v2inv[compmix.pois$igz] <- 1/v2
    yh2[compmix.pois$igz] <- yh[-(1:n)]
    
    dBP <- 1/(diag(prior$invB0) + theta^2*(1/v1 + v2inv))
    bP <- dBP*theta*(yh[1:n]/v1 + yh2*v2inv)
    btilde <- rnorm(max(model$Zp), bP, sqrt(dBP))
    
    # perform a random sign-switch (sign-switching step)
    sswitch <- sign(runif(1) - 0.5)
    theta   <- theta*sswitch
    btilde  <- btilde*sswitch 
  } else {
    btilde <- NULL  
    theta  <- NULL  
  }
  
  ## --- (iii) sample the variance parameter of Student-t/ or normal slab
  ind  <- c(delta, gamma)
  if (model$ri == 1){
    effBeta <- c(t(beta),t(theta))
  } else {
    effBeta <- t(beta)
  }
  
  psi <- draw_psi(effBeta, ind, prior)
  
  # returns updated par-list with parameters used for subsequent step
  return(list(beta = c(mu_beta, beta), delta = delta, pdelta = pdelta, 
              omega = omega, psi = psi, btilde = btilde, theta = theta, 
              gamma = gamma, pgamma = pgamma, pi = pi))
}
