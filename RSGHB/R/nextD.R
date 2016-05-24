nextD <- function(a, b, env)
{		
     `a-b` <- matrix(a, nrow = env$gNP, ncol = env$gNIV, byrow = TRUE) - b
     H <- env$pvMatrix + t(`a-b`) %*% `a-b`
     
     T <- t(chol(solve(H)))
     u <- matrix(rnorm(env$gNIV * (env$gNP + env$gNIV + env$degreesOfFreedom)), nrow = env$gNIV, ncol = env$gNP + env$gNIV + env$degreesOfFreedom)
     S <- (T %*% u) %*% t(T %*% u)
     return(solve(S))
}
