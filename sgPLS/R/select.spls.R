select.spls <- function(model){
  res.select <- list()
  for (h in 1:model$ncomp){
    set.ind.zero <- which(model$loadings$X[,h]!=0)
    names(set.ind.zero) <- model$names$X[set.ind.zero]
    res.select[[h]] <- set.ind.zero
  }
  select.x <- res.select  
  ind.total.x <- sort(unique(unlist(select.x)))
  names(ind.total.x) <- model$names$X[ind.total.x]
  
  res.select <- list()
  for (h in 1:model$ncomp){
    set.ind.zero <- which(model$loadings$Y[,h]!=0)
    names(set.ind.zero) <- model$names$Y[set.ind.zero]
    res.select[[h]] <- set.ind.zero
  }
  select.y <- res.select
  ind.total.y <- sort(unique(unlist(select.y)))
  names(ind.total.y) <- model$names$Y[ind.total.y] 
  return(list(select.X=select.x,select.Y=select.y,select.X.total=ind.total.x,select.Y.total=ind.total.y)) 
}