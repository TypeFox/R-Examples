pairs.ContaminatedMixt <- function(x, criterion = "BIC", ...){
  
  criterion <- match.arg(criterion,.ICnames())
  res <- getBestModel(x,criterion=criterion)$models[[1]]
  
  
  n <- res$n
  groups <- numeric(n)
  for(i in 1:n){    
    if(res$detection[i,2]=="bad")
      groups[i] <- res$G+1
    else
      groups[i] <- res$group[i]    
  }
  pairs(res$X, pch=c(1:res$G, 19)[groups], col=c(2:(res$G+1), 1)[groups], ...)
  
}
