
#integral(from=0, to=1, function(x) a+b^x)
gendilog <- function(a, b, checkparam=TRUE)
{
  if(!(a +1 >0 && b > 0 && a*(1-b) >= 0) && checkparam)
    return(NaN)
  
  if(a == 0)
    return(log(b)/2)
  else if(b == 1)
    return(log(a+1))
  else if(a > 0)
  {
    res <- (dilog(-1/a) - dilog(-b/a))/log(b)+log(a)
  }else
  {
    res <- (dilog(1+b/a) - dilog(1+1/a))/log(b)
    res <- res + log(a+b) - log(-a)*log((a+b)/(a+1))/log(b)
  }
  res
}
