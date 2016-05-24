## somewhat arbitray values
library(Renext)
set.seed(149)

lambda1 <- 1.0
lambda2 <- 1 + 10*rexp(1)

alpha <- 0.3

test <- rep(NA, 100)

for (i in 1:300) {
  
  ## Choose more values in the difficult zones
  p.test <- c(0.1+0.8*sort(runif(50)),
              0.9+0.05*sort(runif(50)),
              0.95 + 0.05*sort(runif(50)))
  
  H.test <- - log(1-p.test)
  
  q.test <- qmixexp2(p.test,
                     prob1 = alpha,
                     rate1 = lambda1,
                     rate2 = lambda2)

  ## the result should be equal to p.test
  
  p2.test <-  pmixexp2(q.test,
                       prob1 = alpha,
                       rate1 = lambda1,
                       rate2 = lambda2)
  
  test[i] <- max(abs(p2.test - p.test))
  
}

stopifnot(max(abs(test)) < 1e-8)
