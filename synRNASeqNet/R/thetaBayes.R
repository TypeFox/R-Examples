thetaBayes <-
function(cellCounts, priorHyperParam = priorHyperParam){
  npprior <- beta_k(priorHyperParam = priorHyperParam,
                            cellCounts = cellCounts)
  B <- npprior$p * npprior$prior
  
  thetak <- (cellCounts + npprior$prior)/(npprior$n + B)
  theta0 <- npprior$prior/(npprior$n + B)
  
  ans <- list(thetak = thetak, theta0 = theta0, p = npprior$p)
  return(ans)
}
