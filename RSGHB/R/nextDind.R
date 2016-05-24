nextDind <- function(a, b, env)
{
     d <- matrix(0, env$gNIV, env$gNIV)
     b <- matrix(t(a), nrow = env$gNP, ncol = env$gNIV, byrow = TRUE) - b
     
     for(k in 1:env$gNIV)
     {
          t <- 1 + t(b[, k]) %*% b[, k]
          s <- sqrt(1/t)[1, 1] * matrix(rnorm(env$gNP + 1), nrow = env$gNP + 1, ncol = 1)
          
          d[k, k] <- solve(t(s) %*% s)
     }
     
     return(d)
}
