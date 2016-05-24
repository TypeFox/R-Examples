select.sgpls <- function(model){
  group.size <- NULL
  result <- list()
  res.select <- list()
  for (h in 1:model$ncomp){
    set.ind.zero <- which(model$loadings$X[,h]!=0)
    names(set.ind.zero) <- model$names$X[set.ind.zero]
    res.select[[h]] <- set.ind.zero
    nb.group <- length(model$ind.block.x)+1
    size.group <- diff(c(0,model$ind.block.x,ncol(model$X)))
    tab.ind <- c(0,model$ind.block.x,ncol(model$X))
    res <- NULL
    for (group in 1:nb.group){
      set.group <- (tab.ind[group]+1):tab.ind[group+1]
      res <- c(res,list(intersect(set.group,set.ind.zero)))
    }
    group.size <- cbind(group.size,unlist(lapply(res,length)))
    result[[h]] <- res
  }
  select.x <- res.select
  group.size <- cbind(size.group,group.size)
  colnames(group.size) <- c("size",paste("comp",1:model$ncomp,sep=""))
  row.names(group.size) <- 1:nb.group
  
  ind.total.x <- sort(unique(unlist(select.x)))
  names(ind.total.x) <- model$names$X[ind.total.x]
  
  group.size.x <- group.size
  result.x <- result 
  
  if(is.null(model$ind.block.y)|is.null(model$keepY)) {
    group.size.y <- NULL
    result.y <- NULL
    select.y <- NULL
    ind.total.y <- NULL
  } else {
    
    group.size <- NULL
    result <- list()
    res.select <- list()
    for (h in 1:model$ncomp){
      set.ind.zero <- which(model$loadings$Y[,h]!=0)
      names(set.ind.zero) <- model$names$Y[set.ind.zero]
      res.select[[h]] <- set.ind.zero
      nb.group <- length(model$ind.block.y)+1
      size.group <- diff(c(0,model$ind.block.y,ncol(model$Y)))
      tab.ind <- c(0,model$ind.block.y,ncol(model$Y))
      res <- NULL
      for (group in 1:nb.group){
        set.group <- (tab.ind[group]+1):tab.ind[group+1]
        res <- c(res,list(intersect(set.group,set.ind.zero)))
      }
      group.size <- cbind(group.size,unlist(lapply(res,length)))
      result[[h]] <- res
    }
    group.size <- cbind(size.group,group.size)
    colnames(group.size) <- c("size",paste("comp",1:model$ncomp,sep=""))
    row.names(group.size) <- 1:nb.group
    
    group.size.y <- group.size
    result.y <- result 
    select.y <- res.select
    ind.total.y <- sort(unique(unlist(select.y)))
    names(ind.total.y) <- model$names$Y[ind.total.y]
  }
  
  
  return(list(group.size.X=group.size.x,select.group.X=result.x,group.size.Y=group.size.y,select.group.Y=result.y,select.X=select.x,select.Y=select.y,select.X.total=ind.total.x,select.Y.total=ind.total.y))
  
}