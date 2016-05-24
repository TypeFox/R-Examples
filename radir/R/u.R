u <-
function(f, x, beta) 
{
  for (i in 1:length(beta)) 
  {
    assign(eval(parse(text = paste0("pars[", i, "]"))), beta[i])
  }
  res <- function(x){eval({x<-x;f})}
  return(res)
}
