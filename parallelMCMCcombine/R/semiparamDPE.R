semiparamDPE<- function( subchain, bandw = rep(1.0,dim(subchain)[1]), anneal = TRUE, shuff = FALSE)
{    
  
  if (length(dim(subchain)) !=3 )
  {
    stop("The subchain must be an array of dimension c(d,sampT,M).")    
  }  
           
  d          <- dim(subchain)[1]  
  sampletotT <- dim(subchain)[2]    
  M          <- dim(subchain)[3]
  
  if (length(bandw) !=d )    
    {
      stop("The bandw must be a vector of length dim(subchain)[1].")    
    } 
      
  if ( M==1 )
  { 
    theta <- array(subchain[,,1], c(d,sampletotT))
    return (theta)
  }
  else
  {      
    # library(mvtnorm)
    
    # permute subposterior samples
    if (shuff == TRUE)
    {
      for (k in 1:M)
        subchain[,,k]<- subchain[,sample(sampletotT),k]
    }  
    
    #  Compute parameters muhatM & sigmahatM, muhatm & sigmahatm
    #  for the semiparametric estimator subposteriors' product 
    #  (Neiswanger section 3.3)
    
    muhatm        <- matrix(NA,d,M)      
    muhatM        <- rep(NA,d)          
    sigmahatm     <- array(NA,c(d,d,M))    
    sigmahatM     <- matrix(NA,d,d)
    sigmahatM.pre <- matrix(NA,d,d)        

    # compute means, covariance and their inverses
    
    if (d==1)  
    {      
      sigmahatm[1,,] <- apply(subchain[1,,], MARGIN=2, FUN = var)
      muhatm[1,]     <- colMeans(subchain[1,,], dims=1)
    }
    else
    {
      for (m in 1:M) { 
        muhatm[,m]     <- rowMeans(subchain[,,m], dims=1)
        sigmahatm[,,m] <- cov(t(subchain[,,m]))
      }
    }
    
    sigmahatm.inverse <- array(NA,c(d,d,M))
    
    # compute the inverses of covariance matrices
    
    for (k in 1:M){
      
      res <- try( sigmahatm.inverse[,,k] <- solve(sigmahatm[,,k]), silent=TRUE)  
      
      if (class(res) == "try-error")
      {  
        stop(paste("Computation of the inverse of a covariance matrix for",
                   "one of the sample vectors in the data subset #",k,"is failed.",
                   "Here is the system R error-message:\n",attr(res,"condition")))  
      }
    }            
    
    sigmahatM.pre <- rowSums(sigmahatm.inverse, dims=2)    
    sigmahatM     <- solve(sigmahatM.pre)
    
    wvec <- rep(0,d)
    
    for (k in 1:M)
    { wvec <- wvec + sigmahatm.inverse[,,k] %*% muhatm[,k] }
    
    muhatM <- sigmahatM %*% wvec            
    
    
    ######## Semiparametric method of Neiswanger p. 6 ######## 
    
    theta        <- matrix(NA,nrow=d,ncol=sampletotT)
    tvector      <- array(NA,M)            
    W_t_denom    <- rep(NA,M)            
    theta_tmat   <- matrix(NA,d,M)
    theta_tmat_m <- array(NA,d)
    theta_tdotbar<- rep(NA,d)    
    theta_cdotbar<- rep(NA,d)      
    sigmatdot    <- matrix(NA,d,d)
    mutdot       <- rep(NA,d)          
    diagh2  <- matrix(NA,d,d)
    diagh2M <- matrix(NA,d,d)
    
    
    # 1. draw {t1,t2,...,t_sampletotT} from unif(1,...,sampletotT)
    #    and pre-compute theta_tmat & theta_tdot         
    
    for (m in 1:M) {
      t_m <- sample(1:sampletotT, size = 1)
      theta_tmat[,m]<- subchain[ ,t_m, m]
      W_t_denom[m]  <- dmvnorm(theta_tmat[,m], mean = muhatm[,m], as.matrix(sigmahatm[,,m]),log=TRUE)
    }
    theta_tdotbar <- rowMeans(theta_tmat, dims=1)
    
    # 2. for i in 1:T do
    
    for (i in 1:sampletotT)
    { 
      
      #3. Set h <- i^(-1/(4+d))
      
      if (anneal == FALSE) { h <- bandw }
      else                 { h <- bandw * i^(-1/(4+d)) }
      
      diagh2  <- diag(h^2,d,d)
      diagh2M <- diag(h^2/M,d,d)
      
      #4. for m=1 to M loop
      
      for (m in 1:M){ 
        
        # 5. set c. <- t. # skip this step because c. and t. will differ only by c_m                    
        
        # 6. draw c_m from unif({1,...,T)}:   
        
        c_m <- sample(1:sampletotT,1)
        
        # 7. draw u from unif[0,1]
        
        u <- runif(1,min=0.0,max=1.0)
        
        # 8. calculate W_c. and W_t., check if u < (W_c./W_t.)          

        # compute W_t.
        
        log_W_t <- sum(dmvnorm(t(theta_tmat), mean = theta_tdotbar, sigma = diagh2,  log=TRUE)) + 
                   dmvnorm(theta_tdotbar, mean = muhatM, sigma = (sigmahatM + diagh2M), log=TRUE) - 
                   W_t_denom[m]
        
        # compute W_c.                
        
        theta_tmat_m   <- theta_tmat[,m] # store temporarily m-th column of theta_t. matrix
                
        theta_tmat[,m] <- subchain[, c_m, m] # theta_t. matrix becomes theta_c. 
        
        theta_cdotbar  <- rowMeans(theta_tmat, dims=1)                     
        
        W_c_denom_m <- dmvnorm(theta_tmat[,m], mean = muhatm[,m], sigma = as.matrix(sigmahatm[,,m]), log=TRUE)
        
        log_W_c <- sum(dmvnorm(t(theta_tmat), mean = theta_cdotbar, sigma = diagh2,  log=TRUE)) + 
                   dmvnorm(theta_cdotbar, mean = muhatM, sigma= sigmahatM + diagh2M, log=TRUE)  - 
                   W_c_denom_m
        
        
        # check if u < ( W_c./ W_t. )
        
        if ( u == 0.0 ) { flag = TRUE }
        else
        {            
          if ( (log_W_t == -Inf) && (log_W_c == -Inf) )
          {
            flag = FALSE
          }
          else
          {            
            if ( log(u) + log_W_t < log_W_c )
            { 
              flag = TRUE
            }
            else { flag = FALSE }          
          }                        
        }
        
        if (flag == TRUE) # if accepted theta_tmat is not changed b/c it has been already modified
        {                     
          theta_tdotbar  <- theta_cdotbar
          
          W_t_denom[m]   <- W_c_denom_m
        } 
        else # if not accepted return the m-th column back
        {
          theta_tmat[,m] <- theta_tmat_m 
        }
        #10. end if 
        
      }  #11. end m
      
      # print(paste("i=",i))
      
      #12. Draw theta_i from  mvnormal(mu_t., sigma_t.)
      
      sigmatdot <- solve( diag(M/h^2,d,d) + sigmahatM.pre )
      
      mutdot    <- sigmatdot %*% ( diag(M/h^2,d,d) %*% theta_tdotbar + sigmahatM.pre %*% muhatM )
      
      theta[,i] <- rmvnorm(n=1, mean=mutdot, sigma=sigmatdot)
      
    }  # 13. end i
    
    return(theta)
  }
  
}