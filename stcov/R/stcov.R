#' @importFrom stats isoreg

sum_diffs <- function(l) {
  rowSums(1 / (outer(l, l, '-') + diag(Inf, length(l))))
}

#' Stein's raw (unisotonized) eigenvalue estimates
#'
#' @param l Sample eigenvalues
#' @param n Number of observations
#' @return Estimated eigenvalues
#' @examples
#' p <- 5
#' n <- 10
#' S <- rWishart(1, n, diag(p))[,,1]
#' l <- eigen(S)$val
#' stein_eig(l, n)
#' @export
stein_eig <- function(l, n) {
  p <- length(l)
  k <- min(n, p)
  l <- l[1:k]
  phi <- l * n / (abs(n - p) + 1 + 2 * l * sum_diffs(l))
  if (k < p) {phi <- c(phi, rep(0, p - k))}
  return(phi)
}

#' Stein's isotonized eigenvalue estimates
#'
#' @param l Sample eigenvalues
#' @param n Number of observations
#' @return Estimated eigenvalues
#' @examples
#' p <- 5
#' n <- 10
#' S <- rWishart(1, n, diag(p))[,,1]
#' l <- eigen(S)$val
#' iso_eig(l, n)
#' @export
iso_eig <- function(l, n) {
  alpha <- l / stein_eig(l, n)

  block <- 1:length(l)

  for (j in length(l):2) {
    if (alpha[j] < 0) {
      l[j - 1] <- l[j - 1] + l[j]
      l <- l[-j]
      alpha[j - 1] <- alpha[j - 1] + alpha[j]
      alpha <- alpha[-j]
      block[j:length(block)] <- block[j:length(block)] - 1
    }
  }

  ratios <- l / alpha
  j <- length(l) - 1
  while (j > 0) {
    if ((length(ratios) >= j + 1) && (ratios[j + 1] >= ratios[j])) {
      l[j + 1] <- l[j + 1] + l[j]
      l <- l[-j]
      alpha[j + 1] <- alpha[j + 1] + alpha[j]
      alpha <- alpha[-j]
      ratios <- ratios[-j]
      ratios[j] <- l[j] / alpha[j]
      block_ind <- which(block > block[j])[1]
      block[block_ind:length(block)] <- block[block_ind:length(block)] - 1
      j <- j + 1
    }
    j <- j - 1
  }

  return(ratios[block])
}

#' Stein/Haff's ordered eigenvalue estimates
#'
#' @param l Sample eigenvalues
#' @param n Number of observations
#' @return Estimated eigenvalues
#' @examples
#' p <- 5
#' n <- 10
#' S <- rWishart(1, n, diag(p))[,,1]
#' l <- eigen(S)$val
#' haff_eig(l, n)
#' @export
haff_eig <- function(l, n) {
  p <- length(l)
  k <- min(n, p)
  phi <- 1 / isoreg(1 / stein_eig(l, n)[1:k])$yf
  if (k < p) {phi <- c(phi, rep(0, p - k))}
  return(phi)
}

#' Stein's isotonized covariance estimator
#'
#' @param S Sample covariance matrix
#' @param n Number of observations
#' @return Estimated covariance matrix
#' @examples
#' p <- 5
#' n <- 10
#' S <- rWishart(1, n, diag(p))[,,1]
#' iso_cov(S, n)
#' @export
iso_cov <- function(S, n) {
  p <- nrow(S)
  decomp <- eigen(S, symmetric=TRUE)
  l <- decomp$val
  H <- decomp$vec
  phi_iso <- iso_eig(l, n)
  return(H %*% (phi_iso * t(H)))
}

#' Stein/Haff's covariance estimator
#'
#' @param S Sample covariance matrix
#' @param n Number of observations
#' @return Estimated covariance matrix
#' @examples
#' p <- 5
#' n <- 10
#' S <- rWishart(1, n, diag(p))[,,1]
#' haff_cov(S, n)
#' @references
#'   Haff, L. R. "The Variational Form of Certain Bayes Estimators." The Annals
#'   of Statistics 19, no. 3 (1991): 1163-1190.
#' @references 
#'   Lin, S.P. and Perlman, M.D.. "A Monte Carlo comparison of four
#'   estimators of a covariance matrix." Multivariate Analysis 6 (1985): 411-429.
#' @references
#'   Stein, C. "Estimation of a covariance matrix". Rietz Lecture (1975).

#' @export
haff_cov <- function(S, n) {
  p <- nrow(S)
  decomp <- eigen(S, symmetric=TRUE)
  l <- decomp$val
  H <- decomp$vec
  phi_haff <- haff_eig(l, n)
  return(H %*% (phi_haff * t(H)))
}
