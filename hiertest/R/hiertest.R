#' Computes contrasts
#'
#' Computes the main effect and interaction contrasts, "w" and "z", respectively.
#'
#' @param x n by p design matrix
#' @param y binary (0 or 1) vector of length n indicating class
#' @param type determines whether Fisher transform should be applied to interaction contrasts.  Default: Fisher
#'
#' @return Returns w and z.
#'
#' @export
compute.contrasts <- function(x, y, type=c("Fisher", "simple")) {
  # compute main effect and interaction contrasts, "w" and "z", respectively.
  n <- nrow(x)
  p <- ncol(x)
  stopifnot(length(y) == n)
  stopifnot(y %in% 0:1)
  stopifnot(table(y) > 3) # need at least 2 observations in each class
  xbar0 <- colMeans(x[y==0, ])
  xbar1 <- colMeans(x[y==1, ])
  s0 <- apply(x[y==0, ], 2, sd)
  s1 <- apply(x[y==1, ], 2, sd)
  n0 <- sum(y==0)
  n1 <- sum(y==1)
  w <- (xbar0-xbar1) / sqrt(s0^2/n0 + s1^2/n1)
  rho0 <- cor(x[y==0, ])
  rho1 <- cor(x[y==1, ])
  if (type[1]=="Fisher") {
    z <- (atanh(rho0) - atanh(rho1)) / sqrt(1/(n0-3) + 1/(n1-3))
    diag(z) <- 0
  }
  else if (type[1]=="simple") {
    xx0 <- scale(x[y==0, ])
    xx1 <- scale(x[y==1, ])
    zz0 <- array(NA, c(n0, p, p))
    zz1 <- array(NA, c(n1, p, p))
    for (j in seq(p))
      for (k in seq(p)) {
        zz0[, j, k] <- xx0[, j] * xx0[, k]
        zz1[, j, k] <- xx1[, j] * xx1[, k]
      }
    m0 <- m1 <- s0 <- s1 <- matrix(NA, p, p)
    for (j in seq(p))
      for (k in seq(p)) {
        m0[j, k] <- mean(zz0[, j, k]) # should be  ((n-1)/n) * rho0
        m1[j, k] <- mean(zz1[, j, k])
        s0[j, k] <- sd(zz0[, j, k])
        s1[j, k] <- sd(zz1[, j, k])
      }
    z <- (m0 - m1) / sqrt(s0^2 / n0 + s1^2 / n1)
  }
  else stop("Not implemented")
  list(w=w, z=z)
}

#' Convex Hierarchical Testing Method
#'
#' This is the main function, implementing the Convex Hierarchical Testing (CHT)
#' procedure.  The CHT procedure produces a set of test statistics for both main
#' effects and interactions with the property that an interaction's statistic is
#' never larger than at least one of its two main effects.  This is accomplished
#' by formulating a convex optimization problem that enforces a hierarchical 
#' sparsity relationship between the main effects and interactions.  The result 
#' is that interactions with large main effects receive a "boost" relative to 
#' those that do not.
#' 
#' The Convex Hierarchical Testing test statistics are the knots of the CHT 
#' optimization problem.  That is, the statistic for a given main effect or 
#' interaction is the value of lambda at which the corresponding parameter 
#' becomes nonzero in the regularization path.  Theorem 1 of the CHT paper gives
#' the closed form expression used to compute these knots (recall that for the
#' interaction test statistics, one takes the maximum of the two corresponding
#' knots).
#' 
#' In Section 2.1 of the CHT paper, the raw main effect and interaction 
#' contrasts are defined.  These are referred to as "w" and "z" in the paper. 
#' The main effect contrast "w" is the standard two-sample t-statistic. The 
#' interaction contrast "z" is the normalized difference of the Fisher 
#' transformed sample correlations between the two classes.  If one instead uses
#' \code{type="simple"}, we simply take for "z" a two-sample statistic on the
#' products of features.  We recommend that \code{type="Fisher"} be used instead
#' of \code{"simple"}.
#'
#' @param x n by p design matrix
#' @param y binary (0 or 1) vector of length n indicating class
#' @param type determines whether Fisher transform should be applied to
#'   interaction contrasts.  See below for explanation. Default is Fisher
#'   and is the recommended choice.
#'   
#' @return A hiertest object, which consists of an ordered list of the
#'   main effects and interactions and a vector indicating which of these are
#'   interactions.
#'   
#' @seealso \code{\link{estimate.fdr}}
#'
#' @references Bien, Simon, and Tibshirani (2015) Convex Hierarchical Testing 
#' of Interactions. Annals of Applied Statistics. Vol. 9, No. 1, 27-42.
#' 
#' @examples
#' 
#' # generate some data accoring to the backward model:
#' set.seed(1)
#' n <- 200
#' p <- 50
#' y <- rep(0:1, each=n/2)
#' x <- matrix(rnorm(n*p), n, p)
#' colnames(x) <- c(letters,LETTERS)[1:p]
#' # make some interactions between several pairs of variables:
#' R <- matrix(0.3, 5, 5)
#' diag(R) <- 1
#' x[y==1, 1:5] <- x[y==1, 1:5] %*% R
#' # and a main effect for variables 1 and 3:
#' x[y==1, 1:5] <- x[y==1, 1:5] + 0.5
#' testobj <- hiertest(x=x, y=y, type="Fisher")
#' # look at test statistics
#' print(testobj)
#' plot(testobj)
#' \dontrun{
#' lamlist <- seq(5, 2, length=100)
#' estfdr <- estimate.fdr(x, y, lamlist, type="Fisher", B=200)
#' plot(estfdr)
#' print(estfdr)
#' # the cutoff lamlist[70] is estimated to have roughly 10% FDR:
#' estfdr$fdr[70]
#' # this allows us to reject this many interactions:
#' nrejected <- estfdr$ncalled[70]
#' # These are the interactions rejected:
#' interactions.above(testobj, lamlist[70])
#' }
#' 
#' @export hiertest
hiertest <- function(x, y, type=c("Fisher", "simple")) {
  cc <- compute.contrasts(x, y, type=type)
  knots <- get.knots(cc$w, cc$z)
  out <- test.statistics(knots$main.knots, knots$int.knots)
  class(out) <- "hiertest"
  out
}

#' Get significant interactions above a threshold
#' 
#'  Returns all interactions whose statistics are above a certain level
#'  
#'  @param testobj output of \link{hiertest}
#'  @param lambda a threshold above which we're interested in interactions
#'  
#'  @export
interactions.above <- function(testobj, lambda) {
  ints <- testobj$stats[testobj$is.int]
  ints[ints > lambda]
}

#' Plot Convex Hierarchical Testing Statistics
#'
#' This plots the main effect and interaction statistics from the CHT procedure.
#' 
#' @param x output of \link{hiertest}.
#' @param ... additional arguments to \code{plot}
#' 
#' @export
plot.hiertest <- function(x, ...) {
  plot(x$stats,col=x$is.int+1, pch=20, 
       main="Test statistics ordered by magnitude",
       ylab="CHT Statistic", ...)
  legend("topright", legend=c("Main effect", "Interaction"), col=1:2, pch=20)
}

#' Print Convex Hierarchical Testing Statistics
#'
#' This prints the top main effect and interaction statistics from the CHT 
#' procedure. Use \link{estimate.fdr} to get an appropriate cutoff.
#' 
#' @param x output of \link{hiertest}
#' @param nshow show the largest \code{nshow} statistics
#' @param ... additional arguments to \code{print.data.frame}
#' 
#' @export
print.hiertest <- function(x, nshow = 20, ...) {
  top <- x$stats[1:nshow]
  print(data.frame(Name = names(top), Statistic = top), row.names = FALSE, ...)
  cat("... List truncated.  Use estimate.fdr to determine a proper cutoff.", 
      fill = TRUE)
}



reject <- function(test.obj, lamlist) {
  # Returns a set of main effects and interactions "rejected"
  # at the provided set of cutoffs.
  #
  # Args:
  #  test.obj: output of test.statistics
  #  lamlist: set of cutoffs
  lamlist <- sort(lamlist) # small to large means we can shrink set...
  nams <- names(test.obj$stats)
  reject <- list()
  for (i in seq(length(lamlist))){
    ii <- test.obj$stats >= lamlist[i]
    test.obj$stats <- test.obj$stats[ii]
    test.obj$is.int <- test.obj$is.int[ii]
    nams <- nams[ii]
    reject[[i]] <- list(int=nams[test.obj$is.int], main=nams[!test.obj$is.int])
  }
  reject
}

get.knots <- function(w, z) {
  # Gets knots using closed form expression derived in paper
  p <- length(w)
  if (nrow(z) != p) browser()
  stopifnot(nrow(z) == p, ncol(z) == p)
  aw <- abs(w)
  diag(z) <- rep(0, p)
  az <- abs(z)
  
  main.knots <- pmax(aw, (aw + apply(az, 1, max)) / 2)
  int.knots <- matrix(NA, p, p)
  for (j in seq(p)) {
    for (k in seq(p)) {
      if (j == k) next
      int.knots[j, k] <- max(aw[j] - sum(pmax(az[j, -j] - az[j, k], 0)), 0) / 2
    }
  }
  int.knots <- pmin(az, az / 2 + int.knots)
  int.knots <- pmax(int.knots, t(int.knots)) # use max(lam_jk, lam_kj)
  
  list(main.knots=main.knots, int.knots=int.knots)
}

test.statistics <- function(main.stats, int.stats) {
  # given a vector of main effect statistics and a matrix of interaction statistics
  # returns a test statistic object.
  #
  # Note: only looks at upper triangle of int.stats
  p <- length(main.stats)
  int.stats <- int.stats[upper.tri(int.stats)]
  o <- order(-c(main.stats,int.stats))
  main.nams <- names(main.stats)
  if (is.null(main.nams)) main.nams <- 1:p
  int.nams <- outer(main.nams, main.nams, paste, sep=":")
  nams <- c(main.nams, int.nams[upper.tri(int.nams)])
  stats <- c(main.stats, int.stats)
  names(stats) <- nams
  stats <- stats[o]
  list(stats=stats, is.int=o>p)
}

#' Estimate FDR
#'
#' Estimates False Discovery Rate (FDR) based on permutation scheme described in
#' CHT paper (reference below).
#' 
#' @param x n by p design matrix
#' @param y binary (0 or 1) vector of length n indicating class
#' @param lamlist a vector of cutoffs for the statistics
#' @param type determines whether Fisher transform should be applied to
#'   interaction contrasts.  Default: Fisher. See \link{hiertest} for more
#'   information.
#' @param B number of permutations
#'   
#' @return A \code{estfdr} object, which consists of \describe{
#'   \item{\code{ncalled}:}{number of interactions called significant at each
#'   cutoff. Set to NA if 0.} \item{\code{null.ncalled}:}{total number, across
#'   all permutations, of (null) interactions rejected at each cutoff}
#'   \item{\code{fdr}:}{estimate of fdr for each cutoff in lamlist. Set to NA if
#'   no interactions are rejected at this cutoff} }
#'   
#' @seealso \code{\link{hiertest}}
#' @references Bien, Simon, and Tibshirani (2015) Convex Hierarchical Testing 
#' of Interactions. Annals of Applied Statistics. Vol. 9, No. 1, 27-42.
#' 
#' @examples
#' 
#' # generate some data accoring to the backward model:
#' set.seed(1)
#' n <- 200
#' p <- 50
#' y <- rep(0:1, each=n/2)
#' x <- matrix(rnorm(n*p), n, p)
#' colnames(x) <- c(letters,LETTERS)[1:p]
#' # make some interactions between several pairs of variables:
#' R <- matrix(0.3, 5, 5)
#' diag(R) <- 1
#' x[y==1, 1:5] <- x[y==1, 1:5] %*% R
#' # and a main effect for variables 1 and 3:
#' x[y==1, 1:5] <- x[y==1, 1:5] + 0.5
#' testobj <- hiertest(x=x, y=y, type="Fisher")
#' # look at test statistics
#' print(testobj)
#' plot(testobj)
#' \dontrun{
#' lamlist <- seq(5, 2, length=100)
#' estfdr <- estimate.fdr(x, y, lamlist, type="Fisher", B=200)
#' plot(estfdr)
#' print(estfdr)
#' # the cutoff lamlist[70] is estimated to have roughly 10% FDR:
#' estfdr$fdr[70]
#' # this allows us to reject this many interactions:
#' nrejected <- estfdr$ncalled[70]
#' # These are the interactions rejected:
#' interactions.above(testobj, lamlist[70])
#' }
#' 
#' @export
estimate.fdr <- function(x, y, lamlist, type=c("Fisher", "simple"), B=100) {
  # estimates the FDR for hiertest
  #
  # Args:
  #  lamlist: a vector of cutoffs for the statistics
  p <- ncol(x)
  nlam <- length(lamlist)
  ut <- upper.tri(diag(p))
  cc <- compute.contrasts(x, y, type=type) # contrasts of actual dataset
  lams <- get.knots(cc$w, cc$z)$int.knots[ut] # hiertest interaction statistics
  # number of interactions rejected at each cutoff:
  ncalled <- rep(NA, nlam)
  for (i in seq(nlam))
    ncalled[i] <- sum(lams >= lamlist[i])
  ncalled[ncalled==0] <- NA # no estimate for FDR when R=0.
  lamstars <- matrix(NA, B, length(lams))
  for (b in seq(B)) {
    ccstar <- compute.contrasts(x, sample(y), type=type) # contrasts of permuted dataset
    # use actual main effect contrasts but permuted interaction contrasts to get
    lamstars[b, ] <- get.knots(cc$w, ccstar$z)$int.knots[ut] # (null) hiertest interaction statistics
  }
  # total number of (null) interactions rejected at each cutoff
  null.ncalled <- rep(NA, nlam)
  for (i in seq(nlam))
      null.ncalled[i] <- sum(lamstars >= lamlist[i])
  out <- list(ncalled=ncalled, null.ncalled=null.ncalled, fdr=pmin(null.ncalled / B / ncalled, 1),
              lamlist=lamlist)
  class(out) <- "estfdr"
  out
}


#' Plot the estimated FDR
#' 
#' @param x output of \link{estimate.fdr}
#' @param ... additional arguments to pass to plot.
#' 
#' @export
plot.estfdr <- function(x, ...) {
  plot(x$ncalled, x$fdr, type="l", ylab="Estimated FDR", 
       xlab="Number of interactions called", ...)  
}

#' Print the estimated FDR
#' 
#' @param x output of \link{estimate.fdr}
#' @param ... additional arguments to pass to print.data.frame.
#' 
#' @export
print.estfdr <- function(x, ...) {
  print(data.frame("Lambda cutoff" = x$lamlist, "Number called" = x$ncalled, "Estimated FDR" = x$fdr), ...)
}
