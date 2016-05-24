# $Id: varset.R,v 1.2 2002/03/26 16:29:15 hothorn Exp $

varset <- function(N, sigma = 0.1, theta = 90, threshold = 0, u = 1:3)
{
  # create U
  U <- matrix(rep(0, 4), ncol = 2)
  U[1, 1] <- u[1]
  U[1, 2] <- u[2]
  U[2, 1] <- u[3]
  U[2, 2] <- (theta-u[1]*u[3])/u[2]
  lambda <- sqrt(U[1, 1]^2 + U[1, 2]^2)
  U[1, ] <- U[1, ]/lambda
  lambda <- sqrt(U[2, 1]^2 + U[2, 2]^2)
  U[2, ] <- U[2, ]/lambda

  e <- matrix(rnorm(2*N, sd = sigma), ncol = 2, byrow = TRUE)
  expl <- matrix(rnorm(2*N), ncol = 2, byrow = TRUE)
  inter <- t(U %*%t(expl) + t(e))
  z <- (inter > threshold)
  resp <- as.factor(ifelse((z[,1] + z[,2]) > 1, 1, 0))
  colnames(expl) <- c("x1", "x2")
  colnames(inter) <- c("y1", "y2")

  result <- list(explanatory = expl, intermediate = inter, response = resp)
  return(result)
}  
