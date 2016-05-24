LogitScale <- function(x, l) {
  #For x = logit(p) and l = log(s), this function returns logit(sp).
  #Return value is -log(exp(-l) + exp(-(x+l)) - 1)
  #which is algebraically equal to logit(sp).
  #Care is taken to maintain precision in the return value.
  
  #Check inputs
  if (!exists("x")) {
    stop("x is not defined")
  }
  if (!is.vector(x)) {
    stop("x is not a vector")
  }
  if (!is.numeric(x)) {
    stop("x is not numeric")
  }
  
  if (!exists("l")) {
    stop("l is not defined")
  }
  if (!is.vector(l)) {
    stop("l is not a vector")
  }
  if (!is.numeric(l)) {
    stop("l is not numeric")
  }
  
  n <- max(length(x), length(l))
  x <- rep(x, length=n)  
  l <- rep(l, length=n)
  
  ## Identify cases with l not large and closer to zero than x+l.
  ## For these set u= -l and v= -(l+1), otherwise reverse the assignment.
  ## If either l or x+l is near zero, -u will be the closest one
  ## and exp(u)-1 will be evaluated as expm1(u) to retain precision.
  ok1 <- l < log(2) & abs(l) < abs(x+l)
  u <- -l - (!ok1)*x
  v <- -l - ok1*x
  ev <- exp(v) # ev is either exp(-(x+l)) or exp(-l) ...
  eumo <- expm1(u) # ... and eumo is the other choice minus one
  
  ## Return value is to be -log(ev+eumo).
  ## Next two lines choose among three ways to calculate it.
  ## First, calculation of l2=log(eumo+ev) anticipates that ev does not dominate
  ## ev+eumo and handles the special case of u large enough to make eumo=Inf.
  l2 <- ifelse(eumo == Inf, pmax(u, v) + log1p(exp(-abs(u - v))), log(eumo + ev))
  
  ## Second, if ev dominates ev+eumo, then result is -v plus a quantity not far from zero.
  out <- -ifelse(v > log(2*abs(eumo)), v + log1p(eumo/ev), l2)
  
  #return result
  return(out)
}
