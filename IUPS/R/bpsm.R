
##########################################################
# Function to implement BPSM
##########################################################

bpsm <- function( Y, t, X, estimand = "ATE", method = "AI", M = 1, L = 1, K = 10000, S = 1000 ) {

  t1 <- proc.time() # Monitor running time
  
  if ( M!= 1 ) { M <- M }  # number of cross group matches requested
  if ( L!= 1 ) { L <- L }  # numebr of within group matches requested
  if ( K!= 10000 ) { K <- K } # number of simulations in BPSM
  if ( S!= 1000 ) { S <- S }  # number of posterior samples in BPSM
  if ( method != "AI" ) { method <- method }  # default is AI. Others include BPSM or Both.
  if (estimand != "ATE") { estimand <- estimand } # default is ATE. Others include ATT for the AI and both ATT and ATC for BPSM.
    
  N <- length(Y)
  D <- ncol(X)
  
  ### Remove missing data:
  N.missing <- 1
  for(i in 1:N){
    for(j in 1:D) { 
      temp <- 1 
      if( is.na(X[i,j]) ) {temp <- 0; break }
    }
    if ( (!is.na(Y[i])) && temp==1 ) N.missing <- c(N.missing, i)
  }
  
  N.missing <- N.missing[-1]
  X <- X[N.missing,]
  t <- t[N.missing]  
  Y <- Y[N.missing]
  N <- length(Y)
  N1 <- sum(t)  # total number of treatment
  N0 <- N - N1  # total number of control

  ### Reshape the data
  dd <- data.frame(t, Y, X)
  dd <- dd[ order(t, decreasing = T) ,] # sort by treatment status
  t <- dd$t  # re-arrange t so that it has the form (1,...,1,0,...,0)
  Y <- dd$Y  # and re-arrange the rest of the data according to t
  Y0 <- Y[t==0]
  Y1 <- Y[t==1]
  X <- NULL; for (i in 1:D) { X <- cbind(X, dd[,i+2] ) }
  Data <- cbind(X, Y) # {pretreatment, outcome} sorted by t
  
  ### Estimate propensity score
  g <- glm( t ~ X - 1, family=binomial(link="logit"))
  p2 <- g$fitted    # estimated propensity scores 
  
  est <- rep(NA, 3)
  se <-  rep(NA, 3)
  
  ### Matching on estimated propensity scores
  meps <- Match( Y=Y, Tr = t, X = p2, estimand = estimand, M = M)  
  ### estimand for the function "Match()": "ATT"(default), "ATE", "ATC"
  est[1] <- meps$est 
  se[1] <-  meps$se  # standard error taking into account matching uncertainty
  var1 <- se[1]^2
  

  ### Direct variance adjustment using AI method
  if(method == "AI" | method == "Both" ) {
    d <- matrix(NA, N, N)
    r <- matrix(NA, N, N)     
    c <- matrix(NA, D, N)
    
    ### ATE
    if(estimand == "ATE"){
      est[2] <- est[1]
      Hset <- sortps(D, N, N1, p2, L, Data )
      c <- matrix(NA, D, N) 
      for (i in 1:N1){
        dat<-Data[Hset[i,2:(2+L)],]
        c[,i] <- cov(dat[,1:D], dat[,D+1]) * (1-p2[i]) /p2[i]   #####here is to calculate c[,i]
      }
      for (i in (N1+1):N){
        dat<-Data[Hset[i,2:(2+L)],]
        c[,i] <- cov(dat[,1:D], dat[,D+1]) * p2[i] /(1-p2[i])   #####here is to calculate c[,i]
      }    
      c <- t(apply(c, 1, mean, na.rm = T))
      VB <- vcov(g)
      var2 <- var1 - c %*% VB %*% t(c)
    }
    
    ### ATT
    if(estimand == "ATT") {
      c1 <- matrix(NA, D, N)
      c2 <- matrix(NA, D, N)
      dtau <- matrix(NA, D, N)
      
      # within group matched IDs
      Hsett <- sortps(D, N, N1, p2, L, Data ) 
      
      # Cross group matched IDs based on ps matching
      matp <- Match( Y=Y, Tr = t, X = p2, estimand = "ATE", M = M)  
      mid_cp <- matrix(matp$index.control, ncol = M, byrow = T)
      mid_tp <- matrix(matp$index.treated, ncol = M, byrow = T)
      
      # Cross group matched IDs based on covariate matching
      mat <- Match( Y=Y, Tr = t, X = X, estimand = "ATE", M = M)  
      mid_c <- matrix(mat$index.control, ncol = M, byrow = T)
      mid_t <- matrix(mat$index.treated, ncol = M, byrow = T)
      
      ### c1: 
      for( i in 1 : N1) { 
        dat <- Data[Hsett[i,2:(2+L)],]      
        c2[,i] <- cov(dat[,1:D], dat[,D+1]) * (1-p2[i])
        
        YHp <- mean( Y[Hsett[i,3:(2+L)]] )
        YJp <- mean(Y[mid_cp[i,]])
        c1[ ,i] <- X[i,]*p2[i]*(1-p2[i])*(YHp - YJp - est[1])
        
        ### dtau over dtheta. Should match on covariates, not PS.
        YJ <- mean(Y[mid_c[i,]])
        dtau[ ,i] <- X[i,] * p2[i]*(1-p2[i]) * (Y[i] - YJ - est[1])
      }
      for( i in (N1+1) : N) {
        dat <-Data[Hsett[i,2:(2+L)],] 
        c2[,i] <- cov(dat[,1:D], dat[,D+1]) * p2[i]^2 / (1-p2[i])
        
        YHp <- mean( Y[Hsett[i,3:(2+L)]] )
        YJp <- mean(Y[mid_tp[i,]])        
        c1[ ,i] <- X[i,]*p2[i]*(1-p2[i])*(YJp - YHp - est[1])
        
        YJ <- mean(Y[mid_t[i,]])
        dtau[ ,i] <- X[i,]*p2[i]*(1-p2[i]) * (YJ - Y[i] - est[1])              
      }
      c1 <- t(apply(c1, 1, mean, na.rm = T)) * N/N1
      c2 <- t(apply(c2, 1, mean, na.rm = T)) * N/N1
      c <- c1 + c2       
      dtau <- t(apply(dtau, 1, mean, na.rm = T)) * N/N1
      
      est[2] <- est[1]
      VB <- vcov(g)
      # var2 <- var1 - c %*% VB %*% t(c) - dtau %*% VB %*% t(dtau)
      ### Check equation (18), there should be be a plus instead of minus before dtau
      var2 <- var1 - c %*% VB %*% t(c) + dtau %*% VB %*% t(dtau)
    }
    
    negV <- 0
    if(var2<0) { se[2] <- NA ; negV <- 1 } else {se[2] <- sqrt(var2)}
    
    if(method == "AI") {
      estimates <- t( rbind(est[c(1,2)], se[c(1,2)]) )
      rownames(estimates) <- c( "Phat", "AI")
      colnames(estimates) <- c( "est", "se" )
      time <- proc.time() - t1
      tab <- list( estimates = estimates, estimand = estimand, time = time)
      return(tab)
    }
 } 

  if(method == "BPSM" | method == "Both" ) {
    ### BPSM
    ET <- c( g$coef, meps$est ) 
    data <- list( "M", "N", "N1", "N0", "D", "Y", "Y0", "Y1", "t", "X", "ET" )
    inits <- function() { list ("bet" = runif(N), "theta" = runif(D), "sigma" = runif(N) ) }
    
    if(estimand=="ATE"){
      parameters <- c( "ate" )
    } 
    
    if(estimand=="ATT"){
      parameters <- c( "att" )
    } 
    
    if(estimand=="ATC"){
      parameters <- c( "atc" )
    } 
    
    bpsm <- jags(data, inits, parameters, model.file = modelpsm, n.iter = K )
    posterior <- bpsm$BUGSoutput$summar # tes: treatment effects
    
    est[3] <- est[1] 
    se[3] <- posterior[1,2]
    
    
    if(method == "BPSM"){
      estimates <- t( rbind(est[c(1,3)], se[c(1,3)]) )
      rownames(estimates) <- c( "Phat", "BPSM" )
      colnames(estimates) <- c( "est", "se" )
      
      time <- proc.time() - t1
      tab <- list( estimates = estimates, estimand = estimand, time = time, sims = K, posterior = posterior )
      return ( tab )
    }
  }
  
  if(method == "Both"){
    estimates <- t( rbind(est, se) )
    rownames(estimates) <- c( "Phat", "AI", "BPSM" )
    colnames(estimates) <- c( "est", "se" )
    
    time <- proc.time() - t1
    tab <- list( estimates = estimates, estimand = estimand, time = time, sims = K, posterior = posterior )
    return ( tab )
  }
} 