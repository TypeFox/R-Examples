if (FALSE) {

library(robustloggamma)
library(boot)

rlg.fun <- function(data, lambda) {  
  temp <- loggammarob(data, control=loggammarob.control(method="oneWL", bootstrap=TRUE, bootstrap.lambda=lambda))
  ## get("i", envir = .GlobalEnv)
  ## i <- i + 1
  ## assign("i", i, envir = .GlobalEnv)
  ## cat("sono qui", i, "\n")
  est <- c(temp$mu, temp$sigma, temp$lambda)
  return(est)
}
     
rlg.rg <- function(data, est) {
  data <- rloggamma(n=length(data), mu=est[1], sigma=est[2], lambda=est[3])
  return(data)
}

set.seed(1234)
x <- rloggamma(20, 0, 1, 0)
temp <- loggammarob(x)
est <- c(temp$mu, temp$sigma, temp$lambda)

i <- 0
system.time(rlg.boot <- boot(x, function(data) rlg.fun(data, est[3]), R = 100, sim = "parametric", ran.gen = rlg.rg, mle = est))

summary(temp)

boot.ci(boot.out = rlg.boot, type = "perc", index = 1)
boot.ci(boot.out = rlg.boot, type = "perc", index = 2)
boot.ci(boot.out = rlg.boot, type = "perc", index = 3)


}
