BayesUFH <-
function(theta, beta, y, X, b, phi, n, betaprior, Sqsigmaprior, betatype, thetatype, Sqsigmatype)
{
     m <- length(theta)
     p <- length(beta)
     if (thetatype == 1)
          theta <- log(theta)
     else
          theta <- log(theta / (1 - theta))
     theta <- c(theta, rnorm((n-1)*m))
     beta <- c(beta, rnorm((n-1)*p))
     s <- rep(0, m)
     if (Sqsigmatype == 0){
          a0 <- Sqsigmaprior[1]
          b0 <- Sqsigmaprior[2]
          Sqsigmav <- 1.0 / rgamma(n, shape = a0 + m / 2, rate = 1)
          result <- .C("BayesUFH", as.double(theta), as.double(beta), as.double(Sqsigmav), 
               as.double(y), as.double(X), as.double(b), as.double(phi), as.integer(n),
               as.integer(m), as.integer(p), as.double(betaprior), as.double(b0), 
               as.integer(betatype), as.integer(thetatype), as.integer(Sqsigmatype),
               as.integer(s))
     }
      else{
          Sqsigmav <- 1.0 / rgamma(n, shape = m / 2 - 1, rate = 1)
          result <- .C("BayesUFH", as.double(theta), as.double(beta), as.double(Sqsigmav), 
               as.double(y), as.double(X), as.double(b), as.double(phi), as.integer(n),
               as.integer(m), as.integer(p), as.double(betaprior), as.double(Sqsigmaprior), 
               as.integer(betatype), as.integer(thetatype), as.integer(Sqsigmatype), 
               as.integer(s))
     }
     if (thetatype == 1)
          result[[1]] <- exp(array(result[[1]], c(m, n)))
     else
          result[[1]] <- 1 / (1 + exp(-array(result[[1]], c(m, n))))
     result[[2]] <- array(result[[2]], c(p, n))
     result[[16]] <- result[[16]] / (n - 1)
     MCMCsample <- list(theta = result[[1]], beta = result[[2]], sigv = result[[3]], 
         theta.rate = result[[16]], type = "UFH")
     MCMCsample
}
