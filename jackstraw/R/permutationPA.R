#' Permutation Parallel Analysis
#'
#' Estimate a number of significant principal components from a permutation test.
#'
#' Adopted from sva::num.sv, and based on Buja and Eyuboglu (1992)
#'
#' @param dat a data matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param B a number (a positive integer) of resampling iterations.
#' @param threshold a numeric value between 0 and 1 to threshold p-values.
#' @param verbose a logical indicator as to whether to print the progress.
#' @param seed a seed for the random number generator.
#'
#' @return \code{permutationPA} returns
#' \item{p}{a list of p-values for significance of principal components}
#' \item{r}{an estimated number of significant principal components based on thresholding p-values at \code{threshold}}
#'
#' @references Buja A and Eyuboglu N. (1992) Remarks on parrallel analysis. Multivariate Behavioral Research, 27(4), 509-540
#' @export permutationPA
permutationPA = function (dat, B = 100, threshold = 0.05, verbose=TRUE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- ncol(dat)
  m <- nrow(dat)

  uu <- fast.svd(dat, tol = 0)
  ndf <- n - 1
  dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
  dstat0 <- matrix(0, nrow = B, ncol = ndf)
  if(verbose==TRUE) message("Estimating a number of significant principal component: ")
  for (i in 1:B) {
      if(verbose==TRUE) cat(paste(i," "))
      dat0 <- t(apply(dat, 1, sample, replace = FALSE))
      uu0 <- fast.svd(dat0, tol = 0)
      dstat0[i, ] <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
  }
  p <- rep(1, n)
  for (i in 1:ndf) {
      p[i] <- mean(dstat0[, i] >= dstat[i])
  }
  for (i in 2:ndf) {
      p[i] <- max(p[(i - 1)], p[i])
  }
  r <- sum(p <= threshold)
  return(list(r = r, p = p))
}
