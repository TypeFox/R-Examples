BayesYC <-
function(theta, beta, y, X, b, phi, n, betaprior, Sqsigmaprior, Recsigmaprior, betatype, Sqsigmatype)
{
     m <- length(theta)
     p <- length(beta)
     theta <- c(theta, rnorm((n-1)*m))
     beta <- c(beta, rnorm((n-1)*p))
     ai <- Recsigmaprior[1:m]
     bi <- Recsigmaprior[(m+1):(2*m)]
     ni <- Recsigmaprior[(2*m+1):(3*m)]
     Sqsigma <- rgamma(n*m, shape = ai + (ni+1) / 2, rate = 1)
     if (Sqsigmatype == 0){
         a0 <- Sqsigmaprior[1]
         b0 <- Sqsigmaprior[2]
         Sqsigmav <- 1.0 / rgamma(n, shape = a0 + m / 2, rate = 1)
         result <- .C("BayesYC", as.double(theta), as.double(beta), as.double(Sqsigmav), 
             as.double(Sqsigma), as.double(y), as.double(X), as.double(b), as.double(phi), as.integer(n),
             as.integer(m), as.integer(p), as.double(betaprior), as.double(b0), as.double(c(bi, ni)), 
             as.integer(betatype), as.integer(Sqsigmatype))
     }
     else{
         Sqsigmav <- 1.0 / rgamma(n, shape = m / 2 - 1, rate = 1)
         result <- .C("BayesYC", as.double(theta), as.double(beta), as.double(Sqsigmav), 
             as.double(Sqsigma), as.double(y), as.double(X), as.double(b), as.double(phi), as.integer(n),
             as.integer(m), as.integer(p), as.double(betaprior), as.double(Sqsigmaprior), as.double(c(bi, ni)), 
             as.integer(betatype), as.integer(Sqsigmatype))
     }
     result[[1]] <- array(result[[1]], c(m, n))
     result[[2]] <- array(result[[2]], c(p, n))
     result[[4]] <- 1.0 / array(result[[4]], c(m, n))
     MCMCsample <- list(theta = result[[1]], beta = result[[2]], sigv = result[[3]], 
         sig2 = result[[4]], type = "YC")
     MCMCsample
}
