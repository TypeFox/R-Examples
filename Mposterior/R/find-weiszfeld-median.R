##' Find M-posterior given list of samples from subset posteriors using Weiszfeld algorithm of Minsker et al. (2014).
##'
##' Samples from subset posteriors are represented as matrices, with
##' the rows representing atoms and the columns indexing
##' dimensions. Length of \code{subsetAtomsList} equals the number of
##' subsets (and posteriors). M-posterior estimates the weights of its
##' atoms given the atoms of subset posteriors using Weiszfeld algorithm
##' as implemented in \code{findWeiszfeldMedian}; \code{sigma} controls the
##' sensitivity of distance between subset posteriors.
##'
##' @title Find M-posterior using samples from subset posteriors
##' @param subsetAtomsList list of samples from subset posteriors; each element of the list must be a matrix
##' @param sigma  parameter of RBF kernel with form exp(- sigma || x - y ||^2)
##' @param maxit  maximum number of iterations estimating M-posterior
##' @param tol    tolerance for declaring convergence of Weiszfeld algorithm
##' @return list containing Weiszfeld wts of subset posteriors, M-posterior atoms, and the estimation set-up
##' @examples
##' set.seed(12345)
##' ## list that contains subset posterior samples from 2-dim Gaussian density
##' subAtomList <- vector("list", 5)
##' subAtomList[[1]] <- cbind(rnorm(100, mean = 1),  rnorm(100, mean = 1))
##' subAtomList[[2]] <- cbind(rnorm(100, mean = -1),  rnorm(100, mean= -1))
##' subAtomList[[3]] <- cbind(rnorm(100, mean = -1),  rnorm(100, mean = 1))
##' subAtomList[[4]] <- cbind(rnorm(100, mean = 1),  rnorm(100, mean = -1))
##' subAtomList[[5]] <- cbind(rnorm(100, mean = 2),  rnorm(100, mean = 2))
##' library(Mposterior)
##' medPosterior <- findWeiszfeldMedian(subAtomList, sigma = 0.1, maxit = 100, tol = 1e-10)
##' medPosterior
##' summary(medPosterior)
##' plot(medPosterior)
##' @author Sanvesh Srivastava \email{sanvesh@@gmail.com}
##' @references
##' Minsker, S., Srivastava, S., Lin, L., and Dunson, D.B. (2014). Robust and Scalable Bayes via a Median of Subset Posterior Measures. \url{http://arxiv.org/abs/1403.2660}
##' @export
findWeiszfeldMedian <- function (subsetAtomsList, sigma = 0.1, maxit = 100, tol = 1e-10) {
  if (length(subsetAtomsList) <= 1 | class(subsetAtomsList) != "list") {
    stop("subset posteriors must be a list of length > 1", call. = FALSE, domain = NULL)
  }
  if (class(subsetAtomsList[[1]]) == "matrix" & length(subsetAtomsList) > 1) {
    medMeasure <- .Call('Mposterior_findWeiszfeldMedian', PACKAGE = 'Mposterior', subsetAtomsList, sigma, maxit, tol)
    medMeasure$sigma <- sigma
    medMeasure$maxit <- maxit
    class(medMeasure) <- c("Mposterior", "list")
  } else {
    stop("subset posteriors must be a list matrices", call. = FALSE, domain = NULL)
  }
  medMeasure
}

##' Print an object of class Mposterior.
##'
##' Print summary of Mposterior class object; ... are ignored.
##' @title Print an object of class Mposterior
##' @param x object of class Mposterior or output of \code{\link{findWeiszfeldMedian}}
##' @param ... extra arguments are ignored currently
##' @return NULL
##' @author Sanvesh Srivastava \email{sanvesh@@gmail.com}
##' @examples
##' set.seed(12345)
##' ## list that contains subset posterior samples from 2-dim Gaussian density
##' subAtomList <- vector("list", 5)
##' subAtomList[[1]] <- cbind(rnorm(100, mean = 1),  rnorm(100, mean = 1))
##' subAtomList[[2]] <- cbind(rnorm(100, mean = -1),  rnorm(100, mean= -1))
##' subAtomList[[3]] <- cbind(rnorm(100, mean = -1),  rnorm(100, mean = 1))
##' subAtomList[[4]] <- cbind(rnorm(100, mean = 1),  rnorm(100, mean = -1))
##' subAtomList[[5]] <- cbind(rnorm(100, mean = 2),  rnorm(100, mean = 2))
##' library(Mposterior)
##' medPosterior <- findWeiszfeldMedian(subAtomList, sigma = 0.1, maxit = 100, tol = 1e-10)
##' print(medPosterior)
##' @export
print.Mposterior <- function (x, ...) {
  natoms <- x$natoms
  df <- data.frame("Subset" = 1:length(natoms), "Natoms" = natoms, "WeiszfeldWts" = round(x$weiszfeldWts, 2))
  cat("Subset posteriors:", "\n")
  print(df)
  cat("M-posterior has", ncol(x$medianAtoms), "dimensions", fill = TRUE)
  cat("RBF kernel with sigma=", x$sigma, fill = TRUE)
  cat("Weiszfeld algorithm converged in", ncol(x$historyWeiszfeldWts), "iterations", fill = TRUE)
  invisible(x)
}

##' Print summary of an object of class Mposterior.
##'
##' Print summary of Mposterior class object; ... are ignored.
##' @title Print an object of class Mposterior
##' @param object object of class Mposterior or output of \code{\link{findWeiszfeldMedian}}
##' @param ... extra arguments are ignored currently
##' @return NULL
##' @examples
##' set.seed(12345)
##' ## list that contains subset posterior samples from 2-dim Gaussian density
##' subAtomList <- vector("list", 5)
##' subAtomList[[1]] <- cbind(rnorm(100, mean = 1),  rnorm(100, mean = 1))
##' subAtomList[[2]] <- cbind(rnorm(100, mean = -1),  rnorm(100, mean= -1))
##' subAtomList[[3]] <- cbind(rnorm(100, mean = -1),  rnorm(100, mean = 1))
##' subAtomList[[4]] <- cbind(rnorm(100, mean = 1),  rnorm(100, mean = -1))
##' subAtomList[[5]] <- cbind(rnorm(100, mean = 2),  rnorm(100, mean = 2))
##' library(Mposterior)
##' medPosterior <- findWeiszfeldMedian(subAtomList, sigma = 0.1, maxit = 100, tol = 1e-10)
##' summary(medPosterior)
##' @author Sanvesh Srivastava \email{sanvesh@@gmail.com}
##' @export
summary.Mposterior <- function (object, ...) {
  natoms <- object$natoms
  df <- data.frame("Subset" = 1:length(natoms), "Natoms" = natoms, "WeiszfeldWts" = round(object$weiszfeldWts, 2))
  cat("Subset posteriors:", "\n")
  print(df)
  cat("M-posterior has", ncol(object$medianAtoms), "dimensions", fill = TRUE)
  cat("RBF kernel with sigma =", object$sigma, fill = TRUE)
  cat("Weiszfeld algorithm converged in", ncol(object$historyWeiszfeldWts), "iterations", fill = TRUE)
  invisible(object)
}

##' Plot an object of class Mposterior.
##'
##' Plot of Mposterior class object; ... are ignored.
##' @param x object of class Mposterior or output of \code{\link{findWeiszfeldMedian}}
##' @param ... extra arguments are ignored currently
##' @return NULL
##' @author Sanvesh Srivastava \email{sanvesh@@gmail.com}
##' @examples
##' set.seed(12345)
##' ## list that contains subset posterior samples from 2-dim Gaussian density
##' subAtomList <- vector("list", 5)
##' subAtomList[[1]] <- cbind(rnorm(100, mean = 1),  rnorm(100, mean = 1))
##' subAtomList[[2]] <- cbind(rnorm(100, mean = -1),  rnorm(100, mean= -1))
##' subAtomList[[3]] <- cbind(rnorm(100, mean = -1),  rnorm(100, mean = 1))
##' subAtomList[[4]] <- cbind(rnorm(100, mean = 1),  rnorm(100, mean = -1))
##' subAtomList[[5]] <- cbind(rnorm(100, mean = 2),  rnorm(100, mean = 2))
##' library(Mposterior)
##' medPosterior <- findWeiszfeldMedian(subAtomList, sigma = 0.1, maxit = 100, tol = 1e-10)
##' plot(medPosterior)
##' @export
plot.Mposterior <- function (x, ...) {
  ww <- x$historyWeiszfeldWts
  nsubs <- nrow(ww)
  plot(ww[1, ], ylim = range(ww), main = "Weiszfeld Wts Updates", xlab = "Iterations", ylab = "Weiszfeld Wts", type = 'l', col = rainbow(nsubs)[1])
  for (ii in 1:nsubs) {
    lines(ww[ii, ], type = 'l', col = rainbow(nsubs)[ii])
  }
  invisible(x)
}
