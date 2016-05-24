change.log <- function(ruleset, node.n, simutry, data)
{ 
  nnn <- ncol(ruleset)  
  loglik.change <- rep(NA, nnn) 
  names(loglik.change) <- colnames(ruleset)
  if (nnn > 1)
  {
    for (j in 1:nnn)
    {
      n.node <- newnode(ruleset, j, data) 
      loglik.change[j] <- delta.log(n.node, node.n, simutry)
    }
    
  }
  else loglik.change<-lik(simutry, data) 
  return(loglik.change)
  
}
