#### Bayesian variable selection for the negative binomial model 
#### using a Dirac spike and Student-t/ or normal slab
#### 
#### (invisible function)
#### last change: 2016/03/24

select_negbin <- function(y, X, offset, mcomp, compmix.pois = NULL, cm1, model, 
                          prior, mcmc, param, imc){
  
  # linear predicor in Poisson model
  muP <- X%*%param$beta
  linp <- muP  
  
  n <- length(y)
  
  #### Step A --- data augmentation for the Poisson model 
  lambda  <- exp(linp)*param$g
  
  
  ## A.1 sample the inter-arrival times of assumed Poisson process in [0,1]
  
  # IAMS - mixture components
  if (model$family == "pogit"){
    compmix.pois <- get_mixcomp(y, mcomp)
    #compmix.pois <- do.call(get_mixcomp, list(y, mcomp))  
  } 
  
  augPois <- iams1_poisson(y, offset*lambda, compmix.pois)
  tau1    <- augPois$t1
  tau2    <- augPois$t2
  
  logMu   <- linp + log(param$g) + log(offset)
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
  gS <- c(param$g, param$g[compmix.pois$igz])
  yS <- (-log(tauS) - mR - log(gS) - log(offsetS))*invSig
  
  Xall <- xS*kronecker(matrix(1, 1, model$d + 1), invSig)

  # inverse prior variance of regression effects (updated)
  invA0 <- diag(c(prior$invM0, 1/param$psi), nrow = model$d + 1)
  
  
  #### Step B --- starts Bayesian variable selection
  if (imc > mcmc$startsel && sum(!model$deltafix) > 0){
    
    ## (1) update mixture weights
    incfix <- sum(param$delta == 1)  
    omega  <- rbeta(1, prior$w['wa0'] + incfix, prior$w['wb0'] + model$d - incfix)
    
    
    ## (2) sample the indicators, the regression coefficients and the scale parameters
    ## --- (i) sample the indicators delta_{beta,j}, gamma_{beta} for the slab component
    
    indic <- draw_indicators_nb(yS, Xall, param$delta, omega, model, prior, invA0)
    delta  <- indic$deltanew 
    pdelta <- indic$pdeltanew
  } else {
    delta  <- param$delta
    pdelta <- param$pdelta
    omega  <- param$omega
  }
  
  ## --- (ii_A) sample the (selected) regression effects 
  index <- c(1, which(delta == 1) + 1) 

  Zsel      <- Xall[, index, drop = FALSE]  # Z*=[1, W^delta,atilde]*sqrt(Sigma^-1)
  dsel      <- length(index)
  invA0_sel <- invA0[index, index, drop = FALSE]
  a0_sel    <- prior$a0[index, , drop = FALSE]
  
  AP    <- solve(invA0_sel + t(Zsel)%*%Zsel)  # A = (A0^-1 + (Z*)'Sigma^-1 Z*)
  aP    <- AP%*%(invA0_sel%*%a0_sel + t(Zsel)%*%yS) # a = A(A0^-1*a0 + (Z*)'Sigma^-1*y
  zetaP <- t(chol(AP))%*%matrix(rnorm(dsel), dsel, 1) + aP 
  
  v1        <- matrix(0, 1, model$d + 1) 
  v1[index] <- t(zetaP)
  mu_beta  <- v1[1]
  
  if (model$d > 0){
    beta  <- v1[2:(model$d + 1)]
  } else {
    beta  <- matrix(0, 1, model$d)   
  }
  
  muP <- X%*%c(mu_beta, beta)
    
  ## --- (iii) sample the variance parameter of Student-t/ or normal slab
  effBeta <- t(beta)
  psi <- draw_psi(effBeta, delta, prior)
  
  
  #### Step C --- 
  ## (a) Sample number of degrees of freedom rho (MH step)
  lambda <- exp(muP)*offset
  rhoMH <- exp(rnorm(1, log(param$rho), prior$eps))
  
  llik1 <- sum(lgamma(y + rhoMH)) - n*lgamma(rhoMH) - sum(lgamma(y + 1)) + n*rhoMH*log(rhoMH) + sum(y*log(lambda)) - sum((y + rhoMH)*log(lambda + rhoMH))
  llik0 <- sum(lgamma(y + param$rho)) - n*lgamma(param$rho) - sum(lgamma(y + 1)) + n*param$rho*log(param$rho) + sum(y*log(lambda)) - sum((y + param$rho)*log(lambda + param$rho))
  
  D1 <- llik1 + dgamma(rhoMH, shape = prior$c0, rate = prior$C0, log = TRUE) 
  D0 <- llik0 + dgamma(param$rho, shape = prior$c0, rate = prior$C0, log = TRUE)
  q1 <- log(rhoMH)
  q0 <- log(param$rho)
  rho.alpha <- min(D1 - q0 - D0 + q1, 0)
  u <- runif(1, 0, 1)
  if(log(u) <= rho.alpha){
    rho <- rhoMH
    rho.acc <- 1
  } else {
    rho <- param$rho
    rho.acc <- 0
  }
    
  ## (b) Sample parameter g
  g <- rgamma(n, shape = (y + rho), rate = (lambda + rho))


  # returns updated par-list with parameters used for subsequent step
  return(list(beta = c(mu_beta, beta), delta = delta, pdelta = pdelta, 
              omega = omega, psi = psi, g = g, rho = rho, rho.acc = rho.acc))
}
