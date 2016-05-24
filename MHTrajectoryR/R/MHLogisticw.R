MHLogisticw <- function(y,
                        xx,
                        w,
                        drugs.names,
                        maxit,
                        alpha,
                        initmodel)
{
  t0 <- proc.time()
  n <- sum(w)
  p <- ncol(xx) 
  
  allbic <- rep(0, maxit)
  
  #Initialisation
  if(missing(initmodel)){
    if (p>200){
      current.model <- rep(0,p)
      current.model[sample(1:p,100)] <- 1
    }else{
      current.model <- rbinom(p, size = 1, prob = rep(1/2, p))
    }
  }else{
    current.model <- initmodel 
  } 
     
  
  
  if(sum(current.model)==0)
    current.fit <- glm(y~1., family = binomial, weights=w)
  else{
    datau <- uniquecombs(cbind(y,xx[,which(current.model==1)]))
    ind <- attr(datau,"index")
    tmp <- cumsum(w[order(ind,decreasing = F)])[cumsum(table(ind))]
    tmp <- tmp - c(0,tmp[-length(tmp)])
    
    current.fit <- glm(datau[,1]~datau[,-1], family = binomial, weights=tmp)
  }
  current.model.bic <-  -current.fit$aic + 2*length(current.fit$coef) - log(n)*length(current.fit$coef)
  
  iter <- 1 
  best.model <- current.model
  best.model.bic <- current.model.bic
  best.fit <- current.fit
  
  allbic[iter] <- current.model.bic
  m.best <- nrow(xx)
  
  while(iter < maxit)
  {
    iter <- iter + 1
    candidate.model <- current.model    
    drug  <- sample(1:p, size = sample(1:alpha, 1, prob=rep(1/alpha, alpha)), prob = rep(1/p,p))
    candidate.model[drug] <- 1 - candidate.model[drug]
    
    #compute candidate.model.bic using the same scheme like the first one
    if(sum(candidate.model)==0){
      fla <- "y~1." 
      candidate.fit <- glm(y~1., family = binomial, weights=w)
    }else{
      datau <- uniquecombs(cbind(y,xx[,which(candidate.model==1)]))
      ind <- attr(datau,"index")
      tmp <- cumsum(w[order(ind,decreasing = F)])[cumsum(table(ind))]
      tmp <- tmp - c(0,tmp[-length(tmp)])
      
      candidate.fit <- glm(datau[,1]~datau[,-1], family = binomial, weights=tmp)
    }
    
    candidate.model.bic <-  -candidate.fit$aic + 2*length(candidate.fit$coef) - log(n)*length(candidate.fit$coef)

    
    #compute of the acceptance probability
    if( (runif(1) < exp(candidate.model.bic - current.model.bic)))
    {
      current.model <- candidate.model
      current.model.bic <- candidate.model.bic 
      current.fit <- candidate.fit
      if(current.model.bic > best.model.bic)
      { 
        best.model <- current.model
        best.model.bic <- current.model.bic
        best.fit <- current.fit
        m.best <- nrow(datau)
      }
    }
    
    allbic[iter] <- current.model.bic
    
  }
  output <- list(all.bic = allbic,
                 best.model = best.model, 
                 best.model.bic = best.model.bic, 
                 best.fit.coef =  best.fit$coef,
                 drugs.names = drugs.names,
                 m.best = m.best,
                 time = (proc.time() - t0)[3])
  
  return(output)
}