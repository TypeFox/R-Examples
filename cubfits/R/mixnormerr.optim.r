### Functions for two components of mixture with measurement errors.

### Negative log likelihood.
mixnormerr.nlogL <- function(theta, X, K, debug = .CO.CT$debug){
  prop <- theta[1:(K - 1)]
  prop <- c(prop, 1 - prop)
  mu <- theta[K:(2 * K - 1)]
  sigma2 <- exp(theta[(2 * K):(3 * K - 1)])
  sigma2.e <- exp(theta[3 * K])
 
  ret <- 0
  for(i.k in 1:K){
    tmp <- dnorm(X, mu[i.k], sqrt(sigma2[i.k] + sigma2.e), log = TRUE) +
           log(prop[i.k])
    ret <- ret + exp(tmp)
  }

  ret <- sum(log(ret), na.rm = TRUE)
  if(debug == 1){
    cat("logL: ", ret, "\n", sep = "")
  } else if(debug > 1){
    cat("logL: ", ret, "\n", sep = "")
    print(theta)
  }

  -ret
} # End of mixnormerr.nlogL().


### Constrained optimization for mixture normal with 2 components. 
optim.mixnormerr.logL <- function(X, PARAM){
  K <- PARAM$K
  theta <- c(PARAM$prop[-K],
             PARAM$mu,
             log(PARAM$sigma2),
             log(PARAM$sigma2.e))

  ### Constrains:
  ###   p_k > 0 for all k
  ###   -p_k > -1 for all k
  ###   p_1 + ... + p_{K - 1} > 0
  ###   -p_1 - ... - p_{K - 1} > -1
  ###   log(sigma_k^2) - log(sigma_r^2) > 0 for all k
  ui <- rbind(cbind(diag(K - 1),
                    matrix(0, nrow = K - 1, ncol = 2 * K + 1)),
              cbind(-diag(K - 1),
                    matrix(0, nrow = K - 1, ncol = 2 * K + 1)),
              c(rep(1, K - 1), rep(0, 2 * K + 1)),
              c(rep(-1, K - 1), rep(0, 2 * K + 1)),
              cbind(matrix(0, nrow = K, ncol = 2 * K - 1),
                    diag(K),
                    rep(-1, K)))

  ci <- c(rep(0, (K - 1)),
          rep(-1, (K - 1)),
          0,
          -1,
          rep(0, K))

  ### Drop duplicates for K = 2.
  if(K == 2){
    ui <- ui[-(3:4),]
    ci <- ci[-(3:4)]
  }

  ### Run constrained optimization.
  ret <- constrOptim(theta, mixnormerr.nlogL, grad = NULL, ui, ci,
                     method = "Nelder-Mead",
                     X = X, K = K)
  ret
} # End of optim.mixnormerr.logL().


### Initial parameters.
init.param <- function(X, K = 2){
  tmp <- rowMeans(X, na.rm = TRUE)
  mu <- as.vector(quantile(tmp, prob = seq(0.20, 0.99, length = K)))
  sigma2 <- rep(var(tmp, na.rm = TRUE) / K, K)
  sigma2.e <- mean(apply(X, 1, var, na.rm = TRUE)) / K
  if(sigma2.e > sigma2[1]){
    sigma2.e <- sigma2[1] / 2
  }

  tmp <- matrix(tmp, nrow = K, ncol = length(tmp), byrow = TRUE)
  prop <- tabulate(apply(abs(tmp - mu), 2, which.min)) / length(tmp)

  PARAM <- list(K = K, prop = prop, mu = mu, sigma2 = sigma2,
                sigma2.e = sigma2.e)
  PARAM
} # End of init.param().


### Convert theta to parameters.
get.param <- function(theta, K = 2){
  prop <- theta[1:(K - 1)]
  prop <- c(prop, 1 - prop)
  mu <- theta[K:(2 * K - 1)]
  sigma2 <- exp(theta[(2 * K):(3 * K - 1)])
  sigma2.e <- exp(theta[3 * K])

  ### Convert to a list.
  PARAM <- list(K = K, prop = prop, mu = mu, sigma2 = sigma2,
                sigma2.e = sigma2.e)
  PARAM
} # End of get.param().


### Main function.
mixnormerr.optim <- function(X, K = 2, param = NULL){
  if(is.null(param)){
    PARAM <- init.param(X, K)
  } else{
    K <- param$K
    PARAM <- param
  }
  tmp <- optim.mixnormerr.logL(X, PARAM)
  PARAM.new <- get.param(tmp$par, K)

  ret <- list(param = PARAM.new,
              param.start = PARAM,
              optim.ret = tmp)
  class(ret) <- "mixnormerr"
  ret
} # End of mixnormerr.optim().


### S3 print method.
my.format <- function(x, digits = max(4, getOption("digits") - 3)){
  paste(formatC(x, format = "f", width = -1, digits = digits), collapse = " ")
} # End of my.format().

print.mixnormerr <- function(x, digits = max(4, getOption("digits") - 3), ...){
  cat("prop = ", my.format(x$param$prop, digits), "\n", sep = "")
  cat("mu = ", my.format(x$param$mu, digits), "\n", sep = "")
  cat("sigma2 = ", my.format(x$param$sigma2, digits), "\n",
      "    sd = ", my.format(sqrt(x$param$sigma2), digits), "\n", sep = "")
  cat("sigma2.e = ", my.format(x$param$sigma2.e, digits), "\n",
      "    sd.e = ", my.format(sqrt(x$param$sigma2.e), digits), "\n", sep = "")
  cat("logL = ", my.format(-x$optim.ret$value, digits),
      ", iter = ", paste(x$optim.ret$counts, collapse = " "),
      ", convergence = ", x$optim.ret$convergence, "\n", sep = "")
  invisible()
} # End of print.mixnormerr().


### For plotting.
dmixnormerr <- function(x, param){
  do.call("c", lapply(x, dmixnormerr.one, param))
} # End of dmixnormerr().

dmixnormerr.one <- function(x, param){
  ret <- 0
  for(i.k in 1:param$K){
    tmp <- dnorm(x, param$mu[i.k], sqrt(param$sigma2[i.k] + param$sigma2.e),
                 log = TRUE) +
           log(param$prop[i.k])
    ret <- ret + exp(tmp)
  }
  ret
} # End of dmixnormerr.one().

