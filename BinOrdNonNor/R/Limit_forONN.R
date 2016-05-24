Limit_forONN <-
function (pvec1, skew1, kurto1) 
{
  validate.plist(list(pvec1), 1)
  
  coef <- Fleishman.coef.NN(skew1, kurto1) 
  
  X <- rnorm(1e+05, 0, 1)
  Z <- rnorm(1e+05, 0, 1)
  ZZ <- cbind(1, Z, Z^2, Z^3)  
  Y <- ZZ%*%as.vector(coef)
  XORD <- ordinalize(pvec1, X)
  max <- cor(XORD[order(XORD)], Y[order(Y)])
  min <- cor(XORD[order(XORD, decreasing = TRUE)], Y[order(Y)])
  rm(X, Y)
  return(c(min, max))
}
