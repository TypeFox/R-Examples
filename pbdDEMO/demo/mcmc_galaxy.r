### Initial library
library(pbdMPI, quietly = TRUE)
comm.set.seed(1234, diff = TRUE)                         # set different stream

### Set configurations. 
I.b <- 1000                                              # burn-in iteration
I.t <- 10                                                # thinning iteration
I.n <- 1000                                              # total samples

### Set data and parameters.
x <- MASS::galaxies                                      # obtain galaxies
n <- length(x)
mu.0 <- mean(x)                                          # prior mean of mu.x
sigma.0 <- sqrt(var(x) / n)                              # prior std of sigma.x
sigma.x <- sqrt(var(x))                                  # data std
x.gbd <- x[get.jid(n)]                                   # devide x

### a(theta.new, theta.org)
acceptance <- function(x.gbd, mu.new, mu.org, sigma = sigma.x){
  pi.org <- dnorm(mu.org, mean = mu.0, sd = sigma.0, log = TRUE)
  q.org <- sum(dnorm(x.gbd, mean = mu.org, sd = sigma, log = TRUE))
  q.org <- allreduce(q.org)                              # synchronize

  pi.new <- dnorm(mu.new, mean = mu.0, sd = sigma.0, log = TRUE)
  q.new <- sum(dnorm(x.gbd, mean = mu.new, sd = sigma, log = TRUE))
  q.new <- allreduce(q.new)                              # synchronize

  min(1, exp(pi.new + q.org - pi.org - q.new))
} # End of acceptance

### Hastings-Metropolis MCMC
ret <- NULL
ret.all <- NULL
mu.org <- rnorm(1, mean = mu.0, sd = sigma.0)
# No need to synchronize if diff = FALSE in comm.set.seed().
mu.org <- bcast(mu.org)
for(i in 1:(I.b + I.t * I.n)){
  mu.new <- rnorm(1, mean = mu.0, sd = sigma.0)
  # No need to synchronize if diff = FALSE in comm.set.seed().
  mu.new <- bcast(mu.new)

  a <- acceptance(x.gbd, mu.new, mu.org) 
  U <- runif(1)
  # No need to synchronize if diff = FALSE in comm.set.seed().
  U <- bcast(U)

  if(U <= a){
    mu.org <- mu.new
  }

  ret.all <- c(ret.all, mu.org)
  if(i > I.b && (i %% I.t == 0)){
    ret <- c(ret, mu.org)
  }
}

if(comm.rank() == 0){
  pdf("galaxy_mcmc.pdf", width = 8, height = 4)
    par(mfrow = c(1, 2))
    plot(ret.all, type = "l", main = "MCMC trace",
         xlab = "mu", ylab = "iterations")
    hist(ret, nclass = 40, main = "posterior of mu", xlab = "mu")
    abline(v = mean(ret), col = 2)
  dev.off()

  cat("Total iterations: ", I.b + I.t * I.n, "\n",
      "Burnin: ", I.b, "\n",
      "Thinning: ", I.t, "\n",
      "Total samples: ", I.n, "\n",
      "Posterior mean: ", mean(ret), "\n",
      "95% credible interval:\n",
      sep = "")
  print(quantile(ret, prob = c(0.025, 0.975)))
}

### Finalize jobs.
finalize()
