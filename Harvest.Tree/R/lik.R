lik <- function(try, data)
{
  lik <- 0 
  for (m in 1:length(try)) 
  {
    lik <- likeli(try[[m]]$total, try[[m]]$active) + lik
  }
  lik<-likeli(nrow(data), sum(data$y)) - lik 
  return(lik)
}
