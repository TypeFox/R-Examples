reduced.eigen <-
function(data)
{
  cor = cor(data)
  diag.1 <- diag(rep(1, ncol(data)))
  min.cor = cor - diag.1
  
  reduced.cor = min.cor + diag(smc(data))
  
  e.values = eigen(reduced.cor)
  
  return(e.values)
}
