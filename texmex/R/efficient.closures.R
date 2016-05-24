# Herewith a bunch of factory functions that are used for making priors


.make.mvn.prior <- function(priorParameters) {
  mean <- priorParameters[[1]]
  cov <- priorParameters[[2]]
  if (!isSymmetric(cov, tol=sqrt(.Machine$double.eps),
                   check.attributes = FALSE)) {
    stop("Covariance must be a symmetric matrix")
  }
  # we should also check for a positive definite matrix
  # but hey
  factors <- svd(as.matrix(cov))
  if (min(factors$d) <= 0) {
    stop("Covariance must be nonsingular.")
  }
  num.vars <- length(mean)
  if (ncol(cov) != num.vars) {
    stop("Covariance and mean must be the same size.")
  }

  logdet <- sum(log(factors$d))
  base.ret <- -(num.vars * log(2 * pi) + logdet) / 2
  stacked <- rbind(t(factors$u), t(factors$v))

  function(param) {
    delta <- param - mean
    prod <- as.numeric(stacked %*% delta)
    first <- prod[1:num.vars]
    last <- prod[(1 + num.vars):(2*num.vars)]
    mahab <- as.numeric(first %*% (last / factors$d))
    base.ret - 0.5 * mahab
  }
}


.make.quadratic.penalty <- function(priorParameters) {
  centre <- priorParameters[[1]]
  cov <- priorParameters[[2]]

  factors <- svd(as.matrix(cov))
  if (min(factors$d) == 0) {
    stop("Singular covariance matrix: quadratic penalty impossible")
  }

  stacked <- rbind(t(factors$u), t(factors$v))
  n.prior <- length(centre)
  function(param) {
    delta <- param - centre
    prod <- as.numeric(stacked %*% delta)
    first <- prod[1:n.prior]
    last <- prod[(1 + n.prior):(2*n.prior)]
    as.numeric(first %*% (last / factors$d))
  }
}

.make.lasso.penalty <- function(priorParameters) {
  centre <- priorParameters[[1]]
  cov <- diag(priorParameters[[2]])
  function(param) {
    sum(abs(param - centre) * cov)
  }
}

.make.dummy.penalty <- function(priorParameters) {
  function(param) {
    0
  }
}



.random.spd.matrix <- function(dimn) {
  # generate a random symmetric positive definite matrix
  rank <- 2 * dimn
  while (TRUE) {
    mat <- matrix(rexp(dimn * rank), nrow=rank)
    cov <- t(mat) %*% mat
    # cov is almost surely (technical sense) SPD
    # but let's check
    res <- eigen(cov, TRUE, only.values=TRUE)
    if (min(res$values) > 0) {
      return(cov)
    }
    # otherwise go round again
  }
}

