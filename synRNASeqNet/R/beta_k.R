beta_k <-
function(priorHyperParam, cellCounts){
  n <- sum(cellCounts)
  p <- range(as.numeric(unlist(dimnames(cellCounts))))
  p <- (p[2] - p[1] + 1)^2
  
  if(priorHyperParam == "Jeffreys") prior <- 1/2 else
    if(priorHyperParam == "BLUnif") prior <- 1 else
      if(priorHyperParam == "Perks") prior <- 1/p else
        if(priorHyperParam == "MiniMax") prior <- sqrt(n)/p else
          stop("Unknown Prior")
  
  ans <- list(prior = prior, n = n, p = p)
  return(ans)
}
