homog.test <- function(resmca,var) {
  vs <- varsup(resmca,var)
  N <- sum(vs$weight)
  a <- sqrt(1/outer(1/vs$weight,1/vs$weight,"+"))
  ncp <- resmca$call$ncp
  dev <- list()
  for(i in 1:ncp) dev[[i]] <- abs(outer(vs$coord[,i],vs$coord[,i],"-"))/sqrt(resmca$eig[[1]][i])
  res <- lapply(dev,function(x) a*sqrt((N-1)/N)*x)
  names(res) <- paste('dim',1:ncp,sep='.')
  return(res)
  }