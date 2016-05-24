
ExhaustiveLogisticw <- function(y,
                                x,
                                w,
                                drugs.names)
{
  t0 <- proc.time()
  n <- sum(w)
  p <- ncol(x) 
  
  allmodel <- matrix(c(0,1), 2, 1)
  if (ncol(x)>1){
    for (j in 2:ncol(x)){
      allmodel <- rbind(cbind(rep(0, nrow(allmodel)), allmodel),cbind(rep(1, nrow(allmodel)), allmodel))
    }
  }
  allbic <- rep(0, nrow(allmodel))
  
  for (k in 1:nrow(allmodel)){
    if (sum(allmodel[k,])>0){
      
      datau <- uniquecombs(cbind(y,x[,which(allmodel[k,]==1)]))
      ind <- attr(datau,"index")
      tmp <- cumsum(w[order(ind,decreasing = F)])[cumsum(table(ind))]
      tmp <- tmp - c(0,tmp[-length(tmp)])      
      current.fit <- glm(datau[,1]~datau[,-1], family = binomial, weights=tmp)
      
    }else{
      current.fit <- glm(y~1., family = binomial, weights=w)
    }
    allbic[k] <-  -current.fit$aic + 2*length(current.fit$coef) - log(n)*length(current.fit$coef)
  }
  
  best.model <- allmodel[which(allbic==max(allbic))[1],]
  m.best <- nrow(uniquecombs(cbind(y,x[,which(best.model==1)])))
  
  if (sum(best.model)>0){
    datau <- uniquecombs(cbind(y,x[,which(best.model==1)]))
    ind <- attr(datau,"index")
    tmp <- cumsum(w[order(ind,decreasing = F)])[cumsum(table(ind))]
    tmp <- tmp - c(0,tmp[-length(tmp)])      
    best.fit <- glm(datau[,1]~datau[,-1], family = binomial, weights=tmp)
    
  }else{
    best.fit <- glm(y~1., family = binomial, weights=w)
  }
  
  output <- list(all.bic = allbic,
                 acceptprob = NULL,
                 accepted.model = NULL,
                 best.model = best.model, 
                 best.model.bic = max(allbic), 
                 best.fit.coef = best.fit$coef,
                 drugs.names = drugs.names,
                 m.best = m.best,
                 time = proc.time() - t0,
                 allmodel=allmodel)
  
  return(output)
}