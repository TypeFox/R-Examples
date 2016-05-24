pgenpois <-
function(q, lambda1, lambda2)
{
  foo <- 0
  for (i in 0:q) 
  {
    foo <- foo + dgenpois(i, lambda1 = lambda1, lambda2 = lambda2)
  }
return(foo)
}
