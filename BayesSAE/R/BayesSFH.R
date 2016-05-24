BayesSFH <-
function(theta, beta, lambda, y, X, phi, li1, li2, num, n, betaprior, Sqsigmaprior, betatype, Sqsigmatype)
{
     m <- length(theta)
     p <- length(beta)
     theta <- c(theta, rnorm((n-1)*m))
     beta <- c(beta, rnorm((n-1)*p))
     lambda <- c(lambda, runif((n-1)))
     li1 <- li1 - 1
     li2 <- li2 - 1
     s <- 0
     if (Sqsigmatype == 0){
          a0 <- Sqsigmaprior[1]
          b0 <- Sqsigmaprior[2]
          Sqsigmav <- 1.0 / rgamma(n, shape = a0 + m / 2, rate = 1)
          result <- .C("BayesSFH", as.double(theta), as.double(beta), as.double(Sqsigmav), 
               as.double(lambda), as.double(y), as.double(X), as.double(phi), as.double(num), as.integer(li1),
               as.integer(li2), as.integer(n), as.integer(m), as.integer(p), as.integer(length(li1)),
               as.double(betaprior), as.double(b0), as.integer(betatype), as.integer(Sqsigmatype),
               as.integer(s))
     }
     else{
          Sqsigmav <- 1.0 / rgamma(n, shape = m / 2 - 1, rate = 1)
          result <- .C("BayesSFH", as.double(theta), as.double(beta), as.double(Sqsigmav), 
               as.double(lambda), as.double(y), as.double(X), as.double(phi), as.double(num), as.integer(li1),
               as.integer(li2), as.integer(n), as.integer(m), as.integer(p), as.integer(length(li1)),
               as.double(betaprior), as.double(Sqsigmaprior), as.integer(betatype), as.integer(Sqsigmatype),
               as.integer(s))
     }
     result[[1]] <- array(result[[1]], c(m, n))
     result[[2]] <- array(result[[2]], c(p, n))
     MCMCsample <- list(theta = result[[1]], beta = result[[2]], sigv = result[[3]],
          lam = result[[4]], lam.rate = result[[19]] / (n-1), type = "SFH")
     MCMCsample
}
