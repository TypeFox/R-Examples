agree <- function(object, givgroup, criterion = "BIC"){

  pmcgd <- object
  
  criterion <- match.arg(criterion,.ICnames())
  best <- getBestModel(object,criterion=criterion)
  res <- best$models[[1]]
  
  
  n <- res$n
  obs <- 1:n
  ind.lab <- res$ind.label
  if(is.null(ind.lab))
    ind.unlab <- obs
  if(!is.null(ind.lab)) 
    ind.unlab <- obs[-ind.lab]
  nunlab <- length(ind.unlab)
  
  #if(length(givgroup) != n) 
  #  stop("'givgroup' must have the same length of the sample on which the PMCGD model was fitted")
  
  groups <- numeric(nunlab)
  cont <- 0
  for(i in ind.unlab){
    cont <- cont+1
    if(res$detection[i,2]=="bad")
      groups[cont] <- "bad points"
    else
      groups[cont] <- paste("",res$group[i],sep=" ")     
  }
  if(is.null(ind.lab))
    return(table(givgroup, groups))
  if(!is.null(ind.lab))
    return(table(givgroup[-ind.lab], groups))
  
}
