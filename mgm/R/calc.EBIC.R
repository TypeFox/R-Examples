

calc.EBIC <- function(fit, y, type, emp_lev, v, n, nadj, gam, X, method, weights) {
  
  # function to calculate deviance
  dev_glmnet <- function(object, ...) 
  {
    dev = object$dev
    nulldev = object$nulldev
    (1 - dev) * nulldev
  }
  
  # calculate LL_Null depending on cat/cont
  vartype <- type[v]
  if(method == 'linear') { vartype <- 'g' }
  
  if(vartype=="g") {
    
    mean.i <- coef(fit, s=1)[1] #mean given by intercept model
    LL_null <- sum(dnorm(y,mean.i,1, log=TRUE)*weights)
    
  } else if(vartype=="e") {
    
    mean.i <- coef(fit, s=1)[1] #mean given by intercept model
    LL_null <- sum(dnorm(y,mean.i,1, log=TRUE)*weights)
    
  } else if(vartype=="p") {
    
    mean.i <- coef(fit, s=1)[1] #mean given by intercept model
    LL_null <- sum(dpois(y,exp(mean.i), log=TRUE)*weights)
    
  } else if(vartype=="c") {
    
    n_cats <- emp_lev[v]
    
    ## Compute LL
    m_respdum <- matrix(NA, n, n_cats) # dummy for data
    m_coefs <- matrix(NA, n, n_cats) # dummy for coefficients
    cats <- unique(y)
    for(catIter in 1:n_cats) { 
      m_respdum[,catIter] <- y==cats[catIter] # dummy matrix for 
      m_coefs[,catIter] <- coef(fit, s=1)[[catIter]][1]
    } 
    
    LL_n <- rep(NA, n)
    for(nIter in 1:n) {
      m_LL_parts <- rep(NA,n_cats+1)
      for(catIter in 1:n_cats) {
        m_LL_parts[catIter] <- m_respdum[nIter, catIter] * m_coefs[nIter, catIter]
      }
      m_LL_parts[n_cats+1] <- log(sum(exp(m_coefs[nIter, ])))
      if(m_LL_parts[n_cats+1]==Inf)  m_LL_parts[n_cats+1] <- 10^5 # avoid Inf value for future calculations
      
      LL_n[nIter] <- sum(m_LL_parts)
    }
    LL_null <- 1/n * sum(LL_n*weights)
  }
  
  # calc LL_sat
  LL_sat <- 1/2 * fit$nulldev + LL_null
  
  # calc LL for all lambdas
  dev <- dev_glmnet(fit) # this is already weighted
  LL <- - 1/2 * dev + LL_sat
  
  n_lambdas <- length(fit$lambda)
  
  # calculation of nonzero neighborhoods
  if(vartype!="c") { #continuous case
    coefs_bin <- as.matrix(coef(fit)[-1,]) != 0 #nonzero?
    n_neighbors <- colSums(coefs_bin)
  }
  if(vartype=="c"){ #categorical case
    m_neighbors <- matrix(0,ncol=n_lambdas, nrow=n_cats)
    coefs_bin <- vector("list", length=n_cats)
    for(ca in 1:n_cats){
      coefs_bin[[ca]] <- as.matrix(coef(fit)[[ca]][-1,]) != 0 #nonzero?
    }
    n_neighbors <- colSums(Reduce('+', coefs_bin)!=0) #rule: a predictor has a nonzero parameter with 1 category of the y, then we have a neighborhood relation
  }
  
  # calc all EBICs
  EBIC_lambda <- -2*LL + n_neighbors * log(nadj) + 2*gam*n_neighbors*log(ncol(X))
  
  EBIC_lambda_nNA <- EBIC_lambda[!is.na(EBIC_lambda)]
  lambda_select <- fit$lambda[!is.na(EBIC_lambda)][which.min(EBIC_lambda_nNA)]
  min_EBIC <- min(EBIC_lambda_nNA)
  
  outlist <- list("EBIC"=min_EBIC, "lambda"=lambda_select)
  
  return(outlist)
  
}


