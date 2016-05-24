### This file contains functions related to ASL or ASL*.
### Reference: Kotz, Kozubowski, and Podgorski (2001).

asl.optim <- function(x){
  init <- asl.optim.init(x)
  if(init[4] == 0){
    warning("Initial values are on the boundary.")
  }

  ### For debugging only.
  # cat("Initial parameters:\n")
  # print(init)
  # init.logL <- asl.logL(init[-3], x)    # drop kappa
  # cat("Initial logL:", -init.logL, "\n", sep = "")

  ### Constrains:
  # theta > min(x)
  # theta < max(x)
  # sigma > 0
  ui <- rbind(c(1, 0, 0),
              c(-1, 0, 0),
              c(0, 0, 1))
  ci <- c(min(x), -max(x), 0)

  ### Optimization.
  theta <- init[-3]    # drop kappa
  tmp <- constrOptim(theta, asl.logL, grad = NULL, ui, ci,
                     method = "Nelder-Mead",
                     x = x)
  theta <- tmp$par[1]
  mu <- tmp$par[2]
  sigma <- tmp$par[3]
  kappa <- (sqrt(2 * sigma^2 + mu^2) - mu) / sqrt(2 * sigma)

  ret <- c(theta, mu, kappa, sigma)
  names(ret) <- c("theta", "mu", "kappa", "sigma")

  ### For debugging only.
  # cat("Final parameters:\n")
  # print(ret)
  # final.logL <- asl.logL(ret[-3], x)    # drop kappa
  # cat("Final logL:", -final.logL, "\n", sep = "")

  ret
} # End of asl.optim().

asl.logL <- function(theta, x){
  -sum(dasl(x, theta = theta[1], mu = theta[2], sigma = theta[3], log = TRUE))
} # End of asl.logL().

asl.optim.init <- function(x){
  ### Implementation of Lemma 3.5.2 (p. 173, KKP 2001.)

  ### Step 1.
  theta <- x
  all.h <- do.call("c", lapply(theta, h.theta, x))
  r <- which.min(all.h)

  ### Setp 2.
  if(r == 1){
    ### theta.hat <= x(1)

    ### This case degenerates to an ASL(theta, mu = mean(x) - x(1), 0) which is
    ### equivalent to a positive exponential distribution.
    ### eqn (3.5.75).
    ret <- c(x[r], x[r] - mean(x), 0, 0)
  } else if(r == length(x)){
    ### theta.hat >= x(n)

    ### This case degenerates to an ASL(theta, mu = x(n) - mean(x), 0) which is
    ### equivalent to a negative exponential distribution.
    ### eqn (3.5.81).
    ret <- c(x[r], mean(x) - x[r], 0, 0)
  } else{
    ### x(1) < theta.hat < x(n)

    ### eqn (3.5.121).
    theta.hat <- x[r]
    alpha.beta <- asl.alpha.beta(theta.hat, x)
    sqrt.alpha.beta <- sqrt(alpha.beta)
    forth.root.alpha.beta <- alpha.beta^{0.25}
    kappa.hat <- forth.root.alpha.beta[2] / forth.root.alpha.beta[1]
    sigma.hat <- sqrt(2) * forth.root.alpha.beta[1] * forth.root.alpha.beta[2] *
                 (sqrt.alpha.beta[1] + sqrt.alpha.beta[2])
    mu.hat <- sigma.hat / sqrt(2) * (1 / kappa.hat - kappa.hat)

    ret <- c(theta.hat, mu.hat, kappa.hat, sigma.hat)
  }

  names(ret) <- c("theta", "mu", "kappa", "sigma")
  ret
} # End of asl.optim.init().

h.theta <- function(theta, x){
  ### eqn (3.5.118).

  alpha.beta <- asl.alpha.beta(theta, x)
  sqrt.alpha.beta <- sqrt(alpha.beta)

  ret <- 2 * log(sqrt.alpha.beta[1] + sqrt.alpha.beta[2]) +
         sqrt.alpha.beta[1] * sqrt.alpha.beta[2]
  ret
} # End of h.theta().

asl.alpha.beta <- function(theta, x){
  ### eqn. (3.5.105).

  tmp <- (x - theta)
  alpha <- sum(tmp[tmp >= 0])
  beta <- -sum(tmp[tmp <= 0])

  ret <- c(alpha, beta) / length(tmp)
  ret
} # End of asl.alpha.beta().

