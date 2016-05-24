Limit_forNN <-
function(skew.vec, kurto.vec) 
{
  coef <- Fleishman.coef.NN(skew.vec, kurto.vec) 
  
  X <- rnorm(1e+05, 0, 1)
  XX <- cbind(1,X,X^2,X^3)
  U <- XX%*%coef[1,]  
  
  Y <- rnorm(1e+05, 0, 1)
  YY <- cbind(1,Y,Y^2,Y^3)
  V <- YY%*%coef[2,]  
  
  max <- cor(U[order(U)], V[order(V)])
  min <- cor(U[order(U, decreasing = TRUE)], V[order(V)])
  rm(U, V)
  return(c(min, max))
  
}
