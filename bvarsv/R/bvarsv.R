##################################################
# VAR regressor matrix
##################################################
makeregs <- cmpfun(function(dat,cut,p){
  
  t <- dim(dat)[1]
  M <- dim(dat)[2]
  K <- M + p*(M^2)
  
  Z <- matrix(NA,t*M,K)
  
  for (i in (cut+p+1):t){
    ztemp <- diag(M)
    for (j in 1:p){        
      xtemp <- t(diag(M) %x% dat[i-j,])
      ztemp <- cbind(ztemp,xtemp)
    }
    Z[((i-1)*M+1):(i*M),] <- ztemp
  }
  
  return(Z[((cut+p)*M+1):(t*M),])
})

##################################################
# Helper function 
##################################################
lastcol <- cmpfun(function(x){
  if (is.matrix(x)){
    out <- x[,dim(x)[2]]
  } else {
    out <- x[length(x)]
  }
  return(out)
})

##################################################
# Another helper function
##################################################
beta.reshape <- cmpfun(function(beta, M, p){
auxa <- matrix(rep(NA, length(beta)), nrow = M)
auxa[1:M] <- beta[1:M]
for (pp in 1:p){
  auxa[,((pp-1)*M+2):(pp*M+1)] <- t(matrix(beta[(M+1+(pp-1)*(M^2)):(pp*(M^2)+M)], nrow = M))
}
return(auxa)
})

##################################################
# Regressor matrix for forecasting
##################################################
makeregs.fc <- cmpfun(function(dat,p){
  
  # Dimensions
  T <- dim(dat)[1]
  M <- dim(dat)[2]
  # Keep only most recent lags which are needed for forecasting
  aux <- rbind(dat[(T-p+1):T,,drop=FALSE],matrix(0,1,M))
  # Make regressors 
  return(makeregs(aux,0,p))
  
})

##################################################
# OLS based prior for BVAR-SV model
##################################################
tsprior <- cmpfun(function(priordat,nlag,ndraws=4000){
  
  M <- dim(priordat)[2]
  K <- M + nlag*(M^2) # nr of beta parameters
  numa <- 0.5*M*(M-1) # nr of VCV elements
  
  yt <- t(priordat[(nlag+1):(dim(priordat)[1]),])  # LHS variable
  tau <- dim(yt)[2] # New sample size
  Zt <- makeregs(priordat,0,nlag)             # Regressor matrix
  
  vbar <- matrix(0,K,K)
  xhy <- matrix(0,K,1)
  
  for (i in 1:tau){
    zhat1 <- Zt[((i-1)*M+1):(i*M),]
    vbar <- vbar + t(zhat1) %*% zhat1
    xhy <- xhy + t(zhat1) %*% matrix(yt[,i],M,1)
  }
  aols <- solve(vbar) %*% xhy
  
  sse2 <- matrix(0,M,M)
  for (i in 1:tau){
    e <- matrix(yt[,i] - Zt[((i-1)*M+1):(i*M),] %*% aols,M,1)
    sse2 <- sse2 + e %*% t(e) 
  }
  hbar <- sse2/tau
  hbarinv <- solve(hbar)
  
  vbar <- matrix(0,K,K)
  for (i in 1:tau){
    zhat1 <- Zt[((i-1)*M+1):(i*M),]
    vbar <- vbar + t(zhat1) %*% hbarinv %*% zhat1
  }
  vbar <- solve(vbar)
  
  achol <- t(chol(hbar))
  ssig <- diag(achol)
  for (i in 1:M){
    for (j in 1:M){
      achol[j,i] <- achol[j,i]/ssig[i]    
    }  
  }
  achol <- solve(achol)
  
  a0 <- matrix(0,numa,1)
  ic <- 1
  for (i in 2:M){
    for (j in 1:(i-1)){
      a0[ic,1] <- achol[i,j]
      ic <- ic +  1    
    }
  }
  
  ssig1 <- log(ssig^2)
  
  hbar1 <- solve(tau*hbar)
  hdraw <- matrix(0,M,M)
  a02mo <- matrix(0,numa,numa)
  a0mean <- matrix(0,numa,1)
  
  for (irep in 1:ndraws){
    hdraw <- solve(rWishart(1,tau,hbar1)[,,1])
    achol <- t(chol(hdraw))
    ssig <- diag(achol)
    for (i in 1:M){
      for (j in 1:M){
        achol[j,i] <- achol[j,i]/ssig[i]
      }  
    }
    achol <- solve(achol)
    
    a0draw <- matrix(0,numa,1)
    ic <- 1
    for (i in 2:M){
      for (j in 1:(i-1)){
        a0draw[ic,1] <- achol[i,j]
        ic <- ic +  1    
      }
    }
    a02mo <- a02mo + a0draw %*% t(a0draw)
    a0mean <- a0mean + a0draw
  }
  
  a02mo <- a02mo / ndraws
  a0mean <- a0mean / ndraws
  a02mo <- a02mo - a0mean %*% t(a0mean)
  
  return(list(B_OLS=aols, VB_OLS = vbar, A_OLS = a0, sigma_OLS = ssig1, VA_OLS = a02mo))
}) 

##################################################
# Simulate VAR(1) with TVP and SV
##################################################
sim.var1.sv.tvp <- cmpfun(function(B0 = NULL, A0 = NULL, Sig0 = NULL, Q = NULL, S = NULL, W = NULL, t = 500, init = 1000){

  if (is.null(B0)) B0 <- cbind(rep(0,2), 0.6*diag(2)) # defaults to bivariate process with zero mean and moderately high persistence
  M <- nrow(B0)
  if (is.null(A0)) A0 <- rep(0, 0.5*M*(M-1))
  if (is.null(Sig0)) Sig0 <- rep(0, M)
  if (is.null(Q)) Q <- 1e-10*as.matrix(diag(length(B0))) # defaults to very small variance (almost no parameter drift)
  if (is.null(S)) S <- 1e-10*as.matrix(diag(length(A0)))
  if (is.null(W)) W <- 1e-10*as.matrix(diag(length(Sig0)))
  
  # Reshape initial parameter vector
  
  aux <- rep(NA,M^2 + M)
  aux[1:M] <- B0[,1]
  for (ll in 1:M){
    aux[(M+1+(ll-1)*M):(M+ll*M)] <- B0[ll,2:(M+1)]    
  }
  B0 <- aux  
  
  # Initialize stuff
  
  npar <- length(B0)
  ystar <- matrix(0,1,M)
  Btsim <- B0
  Btall <- array(NA, c(M,M + 1, t))
  Htall <- array(NA, c(M,M,t))
  corall <- matrix(NA,M*(M-1)*0.5,t)
  sigall <- matrix(NA,M,t)
  Atsim <- A0
  Sigtsim <- Sig0
  Q.c <- t(chol(Q))
  S.c <- t(chol(S))
  W.c <- t(chol(W))
  yall <- matrix(NA,t,M)
  
  for (j in 1:(t+init)){
    
    # Draw Bt
    Btsim <- Btsim + Q.c %*% rnorm(npar)
    
    # Draw At
    Atsim <- Atsim + S.c %*% rnorm(length(Atsim))
    
    # Draw Sigt
    Sigtsim <- Sigtsim + W.c %*% rnorm(M)
    
    # Create the VAR covariance matrix
    capAtsim <- diag(M)
    ic <- 1
    for (k in 2:M){
      capAtsim[k,1:(k-1)] <- t(Atsim[ic:(ic+k-2)])
      ic <- ic + k - 1  
    }
    Hsim <- solve(capAtsim)*diag(as.vector(exp(0.5*Sigtsim)))
    Hsim <- Hsim %*% t(Hsim)
    
    # Draw Yt
    Zsim <- makeregs.fc(ystar,1)
    ysim <-  Zsim %*% Btsim + t(chol(Hsim)) %*% rnorm(M)
    ystar <- rbind(ystar, t(ysim))
    
    # Save post-init draws
    
    if (j > init){
    
    Btall[,,j-init] <- Btsim
    Htall[,,j-init] <- Hsim
    yall[j-init,] <- t(ysim)
    
    }
    
  }

  # Return list
  return(list(data = yall, Beta = Btall, H = Htall))
  
}) 

##############################################################################
# Get parameters of mixture normal approximation to log-chi2 distribution 
##############################################################################
getmix <- cmpfun(function(){
  
  q <- c(0.00730, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.25750)      # probabilities
  m <- c(-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -1.08819) # means
  u2 <- c(5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261)    #variances
  return(list(q=q,m=m,u2=u2))
  
})

##################################################
# Estimate Primiceri BVAR with SV and TVP 
##################################################
bvar.sv.tvp <- cmpfun(function(Y, p = 1, tau = 40, nf = 10, pdrift = TRUE, nrep = 50000, nburn = 5000, thinfac = 10, itprint = 10000, save.parameters = TRUE, k_B = 4, k_A = 4, k_sig = 1, k_Q = 0.01, k_S = 0.1, k_W = 0.01, pQ = NULL, pW = NULL, pS = NULL){

  # Input checks

  if (is.null(ncol(Y))) stop( "Y must be a matrix with at least two columns" )
  if (is.matrix(Y)) if (ncol(Y) < 2) stop( "Y must be a matrix with at least two columns" )
  if (p != as.integer(p)) stop( "p: Number of lags must be in {1, 2, 3,...}" )
  if ((tau + p) >= nrow(Y)) stop( "Too few observations for the given choices of p and tau" )
  if (nf < 1) stop( "Please set nf to 1 or higher ")
  if (nf != as.integer(nf)){
	nf <- round(nf)
	warning( "nf rounded to nearest integer" )
  }
  if ( any(c(k_Q, k_S, k_W) <= 0) | any(c(length(k_Q), length(k_S), length(k_W)) != 1) ) stop( "Prior parameters in k_Q, k_S and k_W must be strictly positive scalars" )
  
  # Default values for prior parameters
  
  if (is.null(pQ)) Qprior <- tau else Qprior <- pQ
  if (is.null(pW)) Wprior <- dim(Y)[2] + 1 else Wprior <- pW
  if (is.null(pS)) Sprior <- (1:(dim(Y)[2])) + 1 else Sprior <- pS
  
  # Matrix dimensions
  
  t <- dim(Y)[1]      # Time series and cross section dimensions
  M <- dim(Y)[2]
  numa <- 0.5*M*(M-1) # Nr of VCV matrix elements
  K <- M + p*(M^2)    # Nr of beta parameters
  
  # Make VAR regressor matrix
  
  Z <- makeregs(Y,tau,p)
  
  # Redefine LHS variables
  
  y <- t(Y[(tau+p+1):t,])
  t <- dim(y)[2]        # new sample size for Bayes analysis
  
  # Numerical parameters for initializing draws 
  
  consQ <- 0.0001
  consS <- 0.0001
  consH <- 0.01
  consW <- 0.0001
  
  # Get OLS inputs for prior parameters
  
  prior <- tsprior(Y[1:(tau+p),],p,2000)                
  
  # Make prior parameters (Gaussian priors)
  
  B_0_prmean <- prior$B_OLS
  B_0_prvar <- k_B*prior$VB_OLS 
  A_0_prmean <- prior$A_OLS
  A_0_prvar <- k_A*prior$VA_OLS 
  
  # log Gaussian prior
  
  sigma_prmean <- prior$sigma_OLS
  sigma_prvar <- k_sig*diag(M)  
  
  # inverse Wishart priors (Q = covariance of B(t), W = covariance of log SIGMA(t), S = covariance of A(t))
  
  Q_prmean <- ((k_Q)^2)*Qprior*prior$VB_OLS  
  Q_prvar <- Qprior
  
  W_prmean <- ((k_W)^2)*Wprior*diag(M)
  W_prvar <- Wprior
  
  S_prmean <- vector("list",M-1)
  S_prvar <- matrix(0,M-1,1)
  ind <- 1
  for (ii in 2:M){
    S_prmean[[ii-1]] <- ((k_S)^2)*Sprior[ii-1]*prior$VA_OLS[((ii-1)+(ii-3)*(ii-2)/2):ind,((ii-1)+(ii-3)*(ii-2)/2):ind];
    S_prvar[ii-1,1] <- Sprior[ii-1]
    ind <- ind + ii  
  }
  
  # Initialize empty matrices for posterior draws
  
  Ht <- matrix(1,t,1) %x% (consH*diag(M))
  Htchol <- matrix(1,t,1) %x% (sqrt(consH)*diag(M))
  Qdraw <- consQ*diag(K)
  Sdraw <- consS*diag(numa)
  Sblockdraw <- vector("list",M-1)
  ijc <- 1
  for (jj in 2:M){
    Sblockdraw[[jj-1]] <- Sdraw[((jj-1)+(jj-3)*(jj-2)/2):ijc,((jj-1)+(jj-3)*(jj-2)/2):ijc]
    ijc = ijc + jj  
  }
  Wdraw <- consW*diag(M)
  Btdraw <- matrix(0,K,t)
  Atdraw <- matrix(0,numa,t)
  Sigtdraw <- matrix(0,M,t)
  sigt <- matrix(1,t,1) %x% (0.01*diag(M))
  statedraw <- 5*matrix(1,t,M)
  Zs <- matrix(1,t,1) %x% diag(M)
  prw <- matrix(0,7,1)
  # new stuff (March 29, 2015)
  nrep2 <- nrep/thinfac # effective number of draws
  # re-design forecast output (save only thinned draws)
  fc.y <- fc.m <- array(0,c(M,nf,nrep2))
  fc.v <- array(0,c(M*(M+1)*0.5,nf,nrep2))  
  if (save.parameters == TRUE){
    # arrays for parameter draws
    Bt.alldraws <- array(0, c(K,t,nrep2))
    Ht.alldraws <- Sigt.alldraws <- At.alldraws <- array(0, c(M,M*t,nrep2))
  } else {
    Bt.alldraws <- Ht.alldraws <- NULL
  }
  
  # Storage matrices for (running) posterior means
  
  Bt_postmean <- matrix(0,K,t)
  At_postmean <- matrix(0,numa,t)
  Sigt_postmean <- matrix(0,M,t)
  Ht_postmean <- matrix(0,t,1) %x% (consH*diag(M))
  Qmean <- matrix(0,K,K)
  Smean <- matrix(0,numa,numa)
  Wmean <- matrix(0,M,M)
  sigmean <- matrix(0,t,M)
  cormean <- matrix(0,t,numa)
  sig2mo <- matrix(0,t,M)
  cor2mo <- matrix(0,t,numa)
  
  # Elements of mixture normal approximation to log chi2
  
  tmp <- getmix()
  q <- tmp$q
  m <- tmp$m
  u2 <- tmp$u2
  
  # Gibbs sampler
  
  prints <- seq(from=nburn,to=(nrep+nburn),by=itprint)
  print(paste(Sys.time(),"-- now starting MCMC"))
  
  # auxiliary counter for saved draws
  aux.ct <- 0
  
  for (irep in 1:(nrep+nburn)){
    
    if (irep %in% prints) print(paste(Sys.time(),"-- now at iteration",irep))
    
    # Draw beta
    
    Btdraw <- carterkohn(y,Z,Ht,Qdraw,K,M,t,B_0_prmean,B_0_prvar)$bdraws  
    
    # Draw Q, the covariance of B(t)
    
    sse_2 <- Btdraw[,2:t] - Btdraw[,1:(t-1)]
    sse_2 <- sse_2 %*% t(sse_2)
    Qdraw <- solve(rWishart(1,t+Q_prvar,solve(sse_2+Q_prmean))[,,1])
    
    # Draw alpha
    
    yhat <- alphahelper(y, Z, Btdraw)
    Zc <- -t(yhat);
    sigma2temp <- t(exp(Sigtdraw))
    
    Atdraw <- {}
    ind <- 1
    for (ii in 2:M){ 
      Atblockdraw <- carterkohn(yhat[ii,,drop=FALSE],Zc[,1:(ii-1),drop=FALSE],sigma2temp[,ii,drop=FALSE],matrix(Sblockdraw[[ii-1]],ii-1,ii-1), ii-1,1,t,A_0_prmean[((ii-1)+(ii-3)*(ii-2)/2):ind,],A_0_prvar[((ii-1)+(ii-3)*(ii-2)/2):ind,((ii-1)+(ii-3)*(ii-2)/2):ind,drop=FALSE])$bdraws
      Atdraw <- rbind(Atdraw,Atblockdraw)    
      ind <- ind + ii
    }
    
    sse_2 <- Atdraw[,2:t] - Atdraw[,1:(t-1)]
    sse_2 <- sse_2 %*% t(sse_2)
    
    ijc <- 1
    for (jj in 2:M){
      Sinv <- solve(sse_2[((jj-1)+(jj-3)*(jj-2)/2):ijc,((jj-1)+(jj-3)*(jj-2)/2):ijc] + S_prmean[[jj-1]])
      Sblockdraw[[jj-1]] <- solve(rWishart(1,t-1+S_prvar[jj-1],Sinv)[,,1])
      ijc <- ijc + jj
    }
    
    # Draw Sigma
    
    capAt <- sigmahelper1(Atdraw, M)
    aux <- sigmahelper2(capAt, yhat, q, m, u2, Sigtdraw, Zs, Wdraw, sigma_prmean, sigma_prvar)
    Sigtdraw <- aux$Sigtdraw
    sigt <- aux$sigt
    
    # Draw W, the covariance of SIGMA(t)
    
    sse_2 <- Sigtdraw[,2:t] - Sigtdraw[,1:(t-1)]
    sse_2 <- sse_2 %*% t(sse_2)
    Wdraw <- solve(rWishart( 1, t + W_prvar, solve(sse_2 + W_prmean))[,,1])
    
    # Create the VAR covariance matrix
    
    aux2 <- sigmahelper3(capAt, sigt)
    Ht <- aux2$Ht
    Htsd <- aux2$Htsd 
    
    
    ##############################################################################################################
    # Save post-burnin draws and forecasts
    ##############################################################################################################
    
    if ( (irep > nburn) & ( (irep-nburn) %% thinfac == 0) ){
	  aux.ct <- aux.ct + 1 # counter for saved draws
      
      Bt_postmean <- Bt_postmean + Btdraw 
      At_postmean <- At_postmean + Atdraw 
      Sigt_postmean <- Sigt_postmean + Sigtdraw
	  Ht_postmean <- Ht_postmean + Ht
      Qmean <- Qmean + Qdraw 
      ikc <- 1
      for (kk in 2:M){
        Sdraw[((kk-1)+(kk-3)*(kk-2)/2):ikc,((kk-1)+(kk-3)*(kk-2)/2):ikc] <- Sblockdraw[[kk-1]]
        ikc <- ikc + kk    
      }
      Smean <- Smean + Sdraw 
      Wmean <- Wmean + Wdraw
      
      # Get time-varying correlations and variances
      
      aux3 <- getvc(Ht)
      sigmean <- sigmean + aux3$out1
      cormean <- cormean + aux3$out2    
      
      # Forecasting part
      
	  if (pdrift == TRUE){
	  
      # Forecasts (with parameter drift)
      
      tempfc <- getfcsts(lastcol(Btdraw), lastcol(Atdraw), lastcol(Sigtdraw), Qdraw, Sdraw, Wdraw, t(y), nf, p)   
	  fc.m[,,aux.ct] <- tempfc$mean
	  fc.v[,,aux.ct] <- tempfc$variance
	  fc.y[,,aux.ct] <- tempfc$draw
      
	  } else {
	  
      # Forecasts (without parameter drift)
      
      parmat <- beta.reshape(Btdraw[,t], M, p)       # Current Beta parameter (in matrix format)
      fcdat <- matrix(t(y[,(t-p+1):t]), nrow = p)   # Conditioning data
      auxb <- solve(sigmahelper1(matrix(lastcol(Atdraw), ncol = 1), M)) * diag(as.vector(exp(0.5*lastcol(Sigtdraw))))
      tmpvar <- auxb %*% t(auxb) # Forecast variance

      for (hhh in 1:nf){
        auxl <- varfcst(parmat, tmpvar, fcdat, hhh)
        fc.m[,hhh,aux.ct] <- auxl$mean
        fc.v[,hhh,aux.ct] <- vechC(auxl$variance)
		fc.y[,hhh,aux.ct] <- mvndrawC(auxl$mean, auxl$variance)
      }      
	  
	  }
	  
	  # Save parameter draws
	
	  if (save.parameters == TRUE){
	
	    Bt.alldraws[,,aux.ct] <- Btdraw
		Ht.alldraws[,,aux.ct] <- t(Ht)
		Sigt.alldraws[,,aux.ct] <- Sigtdraw
		At.alldraws[,,aux.ct] <- t(capAt)
	  
	  }
      
    }
    
  }
  
  # Final steps
  
  Bt_postmean <- Bt_postmean / nrep2
  At_postmean <- At_postmean / nrep2
  Sigt_postmean <- Sigt_postmean / nrep2
  Ht_postmean <- Ht_postmean / nrep2
  Qmean <- Qmean / nrep2
  Smean <- Smean / nrep2
  Wmean <- Wmean / nrep2
  sigmean <- sigmean / nrep2
  cormean <- cormean / nrep2
  
  # Prepare outputs
  
  beta.out <- array(NA,c(M,M*p+1,t))
  h.out <- array(NA,c(M,M,t))
  
  for (jj in 1:t){
	beta.out[,,jj] <- beta.reshape(Bt_postmean[,jj], M, p)
	h.out[,,jj] <- Ht_postmean[((jj-1)*M+1):(jj*M),]
  }  
  
  # Return list of outputs
  
  return(list(Beta.postmean = beta.out, H.postmean = h.out, Q.postmean = Qmean, S.postmean = Smean, W.postmean = Wmean, fc.mdraws = fc.m, fc.vdraws = fc.v, fc.ydraws = fc.y, Beta.draws = Bt.alldraws, 
              H.draws = Ht.alldraws, logs2.draws = Sigt.alldraws, A.draws = At.alldraws, M = M, p = p))
  
})

# Helper function to revert vech operator (needed to spot variance elements from VCV matrix)
vels <- cmpfun(function(n){
  aux <- matrix(0, n, n)
  aux[lower.tri(aux, diag = TRUE)] <- 1:(0.5*n*(n+1))
  return(diag(aux)) 
})

predictive.density <- cmpfun(function(fit, v = 1, h = 1, cdf = FALSE){
  # Retrieve number of variables/horizons from fit
  nv <- length(fit$fc.ydraws[, 1, 1])
  nh <- length(fit$fc.ydraws[1, , 1])
  # Input check
  if (v > nv | h > nh | v != round(v) | h != round(h))  stop("Please choose appropriate variable and horizon indices")
  # Locate variance 
  v.ind <- vels(nv)[v]
  # Get means and standard deviations
  m <- fit$fc.mdraws[v, h, ]
  s <- sqrt(fit$fc.vdraws[v.ind, h, ])
  # Finally: Return forecast pdf/cdf
  if (cdf == FALSE){
    return(function(q) sapply(q, function(o) mean(dnorm(o, mean = m, sd = s))))
  } else {
    return(function(q) sapply(q, function(o) mean(pnorm(o, mean = m, sd = s))))
  }
})

predictive.draws <- cmpfun(function(fit, v = 1, h = 1){
  # Retrieve number of variables/horizons from fit
  nv <- length(fit$fc.ydraws[, 1, 1])
  nh <- length(fit$fc.ydraws[1, , 1])
  # Input check
  if (v > nv | h > nh | v != round(v) | h != round(h))  stop("Please choose appropriate variable and horizon indices")
  # Draws for y
  y <- fit$fc.ydraws[v, h, ]
  # Draws for m
  m <- fit$fc.mdraws[v, h, ]
  # Draws for v
  v <- fit$fc.vdraws[vels(nv)[v], h, ] 
  # Return draws in list
  return(list(y = y, m = m, v = v))  
})

parameter.draws <- cmpfun(function(fit, type = "lag1", row = 1, col = 1){
  # Input checks
  if (is.null(fit$Beta.draws)) stop("Parameter draws in fit have not been saved -- please run again, using the option `save.parameters = TRUE'.")
  nms <- c("intercept", paste0("lag", 1:fit$p), "vcv")
  if (!type %in% nms) stop(paste("Type unknown -- please enter one of the following:", paste(nms, collapse = "; ")))
  if (type == "intercept"){
    if (!row %in% 1:fit$M) stop(paste("Invalid row selected -- please enter one of the following:", paste(1:fit$M, collapse = "; ")))
  } else if (type %in% c(paste0("lag", 1:fit$p), "vcv")){
	if (!row %in% 1:fit$M | !col %in% 1:fit$M) stop(paste("Invalid row or col selected -- please enter one of the following:", paste(1:fit$M, collapse = "; ")))
  }
  
  # Intermediate stuff
  Beta.draws <- fit$Beta.draws  
  t <- dim(Beta.draws)[2]
  
  # Get parameters 
  if (type == "vcv"){
    aux.seq <- seq(from = col, by = fit$M, length.out = t)
    out <- t(fit$H.draws[row, aux.seq, ])
  } else {
    if (type == "intercept"){
	  # Parameter output
	  out <- t(Beta.draws[row, , ])
	  # Needed for checking below
      aux.col <- 1
	} else if (type %in% paste0("lag", 1:fit$p)){
	  # Lag order
	  lo <- as.numeric(substr(type, 4, nchar(type)))
	  # Find position in storage output
	  ind <- fit$M + (lo - 1)*(fit$M^2) + (row - 1)*fit$M + col
	  # Get draws
	  out <- t(Beta.draws[ind, , ])
	  # Needed for checking below
	  aux.col <- 1 + (lo-1)*fit$M + col
	}
    # Check (compare to slower beta.reshape function, for one randomly chosen MCMC draw and time period)
	aux.ind <- sample.int(dim(Beta.draws)[3], 1)
	aux.t <- sample.int(dim(Beta.draws)[2], 1)
	if (out[aux.ind, aux.t] != beta.reshape(Beta.draws[, aux.t, aux.ind], fit$M, fit$p)[row, aux.col]) stop("Wrong output...") 
  }
  # Message
  if (type == "intercept"){
    message(paste("Element", row, "of intercept vector. Format: MCMC draws in rows, time in columns."))
  } else {
    message(paste0("Element [", row, ",", col, "] of ", type, " matrix. Format: MCMC draws in rows, time in columns."))
  }
  return(out)
})



ARtoMA <- cmpfun(function(A, nhor){
  # Infer dimensions from A parameters
  M <- nrow(A)
  p <- ncol(A)/M
  Phi <- matrix(0, M, M*nhor) 
  # See recursive formula for MA (see e.g. page 28 in http://www.jmulti.de/download/help/var.pdf)
  for (s in 1:nhor){
    tmp <- matrix(0, M, M)
	for (j in 1:s){
	  if (s == j){
	    aux <- diag(M)
	  } else {
	    aux <- Phi[, ((s-j-1)*M+1):((s-j)*M)]
	  }
	  if (j <= p) tmp <- tmp + aux %*% A[,((j-1)*M+1):(j*M)]
	}
	Phi[,((s-1)*M+1):(s*M)] <- tmp
  }  
  Phi  
})

matmult <- cmpfun(function(mat, mult){
  if (mult == 0){
    out <- diag(nrow(mat))
  } else if (mult == 1){
    out <- mat
  } else {
    out <- mat
    for (c in 1:(mult-1)) out <- out %*% mat    
  }
  return(out)
})

IRFmats <- cmpfun(function(A, H.chol = NULL, nhor, orthogonal = TRUE){
  # Dimensions
  M <- nrow(A) # nr of variables
  
  # AR matrices
  p <- ncol(A)/M
      
  # Construct companion form matrices
  
  if (p > 1){
      
    # VAR matrices
    Ac <- matrix(0, M*p, M*p)
    Ac[1:M, ] <- A
    Ac[-(1:M), 1:(M*(p-1))] <- diag(M*(p-1))
    
  } else {
  
    Ac <- A
  
  }
  
  # Matrix with impulse responses
  Phi <- matrix(0, M, M*(nhor+1)) 
  for (s in 0:nhor){
    aux <- matmult(Ac, s)[1:M, 1:M] 
	if (!is.null(H.chol)) aux <- aux %*% H.chol
	Phi[,(s*M+1):((s+1)*M)] <- aux
  } 
  Phi
  
})

impulse.responses <- cmpfun(function(fit, impulse.variable = 1, response.variable = 2, t = NULL, nhor = 20, scenario = 2, draw.plot = TRUE){

  # Get coefficient draws from fit
  Beta.draws <- fit$Beta.draws
  nd <- dim(Beta.draws)[3]
  if (is.null(t)) t <- dim(Beta.draws)[2]
  out <- matrix(0, nd, nhor + 1)
  M <- fit$M
  p <- fit$p
  
  # Get relevant stuff for VCV matrix
  H.sel <- fit$H.draws[,((t-1)*M+1):(t*M),]
  A.sel <- fit$A.draws[,((t-1)*M+1):(t*M),]
  
  if (scenario == 3){
    # Compute mean of sigma (over both time and MCMC draws, as in Primiceri's code)
    sig <- apply(exp(0.5*fit$logs2.draws), 1, mean)
    sig <- diag(sig)   
  } else {
    sig <- NULL
  }
    
  # Compute IR for each MC iteration
  for (j in 1:nd){
  
    # Cholesky of VCV matrix, depending on specification
    if (scenario == 1){ # No orthogonalization at all
	  H.chol <- NULL
	} else if (scenario == 2){ # Standard orthogonalization
	  H.chol <- t(chol(H.sel[,,j]))
	} else if (scenario == 3){ # Orthogonalization as in Primiceri
	  H.chol <- t(solve(A.sel[,,j])) %*% sig
	}

	# Compute Impulse Responses
    aux <- IRFmats(A = beta.reshape(Beta.draws[,t,j], M, p)[,-1], H.chol = H.chol, nhor = nhor)
	aux <- aux[response.variable, seq(from = impulse.variable, by = M, length = nhor + 1)]  	
	
  	out[j,] <- aux
  }
  # Make plot
  if (draw.plot){
    pdat <- t(apply(out[,-1], 2, function(z) quantile(z, c(0.05, 0.25, 0.5, 0.75, 0.95))))
    xax <- 1:nhor
	matplot(x = xax, y = pdat, type = "n", ylab = "", xlab = "Horizon", bty = "n", xlim = c(1, nhor))
	polygon(c(xax, rev(xax)), c(pdat[,5], rev(pdat[,4])), col = "grey60", border = NA)
	polygon(c(xax, rev(xax)), c(pdat[,4], rev(pdat[,3])), col = "grey30", border = NA)
	polygon(c(xax, rev(xax)), c(pdat[,3], rev(pdat[,2])), col = "grey30", border = NA)
	polygon(c(xax, rev(xax)), c(pdat[,2], rev(pdat[,1])), col = "grey60", border = NA)
	lines(x = xax, y = pdat[,3], type = "l", col = 1, lwd = 2.5)
	abline(h = 0, lty = 2)
  }
  list(contemporaneous = out[,1], irf = out[,-1])  
})
  