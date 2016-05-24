#' @name COV
#' @title \bold{cov}ariance matrix for survival data
#'
#' @include ten.R
#' @include print.R
#' @include asWide.R
#' 
#' @rdname COV
#' @export
#'
COV <- function(x, ...) UseMethod("COV")
#'
#'
#' @param x A \code{numeric} vector of
#'  \emph{number of events}, \eqn{e_t}{e[t]}.
#'  These are assumed to be ordered by discrete times.
#'  \cr
#'  A method is available for objects of \code{class} \code{ten}.
#' @param ... Additional arguments (not implemented).
#' @param reCalc Recalcuate the values?
#'  \cr
#' If \code{reCalc=FALSE} (the default) and the \code{ten} object already has
#' the calculated values stored as an \code{attribute},
#' the value of the \code{attribute} is returned directly.
#'  \cr \cr
#' \bold{--Arguments for the numeric method:}
#' @param n \bold{n}umber at risk (total).
#' @param ncg \bold{n}umber at risk, per \bold{c}ovariate \bold{g}roup.
#'  \cr
#' If there are \eqn{2} groups, this can be given as a \code{vector} with
#' the number at risk for group \eqn{1}.
#' \cr
#' If there are \eqn{\geq 2}{>= 2} groups, it is
#' a \code{matrix} with one column for each group.
#'
#' @details Gives variance-covariance matrix for comparing survival
#' data for two or more groups.
#' \cr
#' Inputs are vectors corresponding to observations at a set of discrete
#' time points for right censored data, except for \eqn{n1},
#' the no. at risk by predictor.
#' \cr
#' This should be specified as a vector for one group,
#' otherwise as a matrix with each column corresponding to a group.
#'
#' @return An \code{array}.
#' \cr
#' The first two dimensions = the number of covariate groups \eqn{K},
#' \eqn{k = 1, 2, \ldots K}.
#' This is the square matrix below.
#' \cr
#' The third dimension is the number of observations
#' (discrete time points).
#' \cr \cr
#' To calculate this, we use \code{x} (= \eqn{e_t}{e[t]} below) and
#' \eqn{n_1}{n1}, the number at risk in covariate group \eqn{1}.
#' \cr
#' Where there are \eqn{2} groups, the resulting sparse square matrix
#' (i.e. the non-diagonal elements are \eqn{0})
#' at time \eqn{t} has diagonal elements:
#'  \deqn{cov_t = - \frac{n_{0t} n_{1t} e_t (n_t - e_t)}{n_t^2(n_t-1)}}{
#'        cov[t] = - n0[t] * n1[t] * e[t] * (n[t] - e[t]) /
#'                  (n[t]^2 * (n[t] - 1))}
#' For \eqn{\geq 2}{>=2} groups, the resulting square matrix
#' has diagonal elements given by:
#'  \deqn{cov_{kkt} = \frac{n_{kt}(n_t - n_{kt}) e_t(n_t - e_t)}{
#'                          n_t^2(n_t - 1)}}{
#'    cov[k, k, t] = n[k, t] * (n[t] - n[k, t]) * e[t] * (n[t] - e[t]) /
#'                   (n[t]^2 * (n[t] - 1))}
#' The off diagonal elements are:
#' \deqn{cov_{klt} = \frac{-n_{kt} n_{lt} e_t (n_t-e_t) }{
#'                         n_t^2(n_t-1)}}{
#'       cov[k, l, t] = - n[k, t] * n[l, t] * e[t] * (n[t] - e[t]) /
#'                      n[t]^2 * (n[t] - 1)}
#' 
#' @note Where the is just one subject at risk \eqn{n=1} at 
#' the final timepoint, the equations above may produce \code{NaN}
#' due to division by zero. This is converted to \code{0} for 
#' simplicity.
#'
#' @seealso Called by \code{\link{comp}}
#' @seealso The name of the function is capitalized 
#' to distinguish it from:
#'  \cr
#' ?stats::cov
#'
#' @keywords survival
#'
#' @rdname COV
#' @method COV ten
#' @aliases COV.ten
#' @export
#' @examples
#' ## Two covariate groups
#' ## K&M. Example 7.2, pg 210, table 7.2 (last column).
#' data("kidney", package="KMsurv")
#' k1 <- with(kidney,
#'            ten(Surv(time=time, event=delta) ~ type))
#' COV(k1)[COV(k1) > 0]
#' ## Four covariate groups
#' ## K&M. Example 7.6, pg 217.
#' data("larynx", package="KMsurv")
#' l1 <- ten(Surv(time, delta) ~ stage, data=larynx)
#' rowSums(COV(l1), dims=2)
#' 
COV.ten <- function(x, ..., reCalc=FALSE) {
    if (!reCalc & !is.null(attr(x, "COV"))) return (attr(x, "COV"))
    ## no. of groups
    g1 <- attr(x, "ncg")
    if (g1 <= 1) stop ("Only valid if more than one covariate group")
    ## if 2 groups only
    if (g1==2) {
        res2 <- x[, (ncg / n) * (1 - (ncg / n)) * ((n - e) / (n - 1)) * e, by=list(t, cg)]
        res2 <- data.table::setkey(res2[, sum(V1), by=t], t)
        res1 <- res2[, V1]
        if (is.nan(res1[length(res1)])) res1[length(res1)] <- 0
        names(res1) <- res2[, t]
    }
    if (g1 > 2) {
        ## same as used in asWide.R
        ## get no. at risk for each unique time and covariate group
        t1 <- data.table::data.table("t" = x[, sort(unique(t))])
        cg1 <- seq.int(attr(x, "ncg"))
        ## abbreviate function
        abbFn <- if (attr(x, "abbNames")) identity else as.integer
        n1 <- sapply(cg1, FUN=function(cg1){
            r1 <- data.table::setkey(x[abbFn(cg)==cg1, ncg, by=t], t)
            r1 <- r1[t1, roll=-Inf]
            data.table::set(r1, i=which(is.na(r1$ncg)), j="ncg", value=0)
            r1[, ncg]
        })
        ## total no. events, no. at risk at each time
        x1 <- x[, list("e"=sum(e), "n"=max(n)), by=t]
        data.table::setkey(x1, t)
        ## 'base variance'; used in all calcuations below
        bv1 <- x1[, e * (n - e) / (n^2 * (n - 1))]
        ## diagonal elements
        r1 <- bv1 * t(apply(n1, MARGIN=1, FUN=
                            function(i) (i * (sum(i) - i))))
        ## off-diagonal elements
        r2 <- bv1 * - t(apply(n1, MARGIN=1, FUN=
                              function(i) apply(utils::combn(i, 2L), MARGIN=2, FUN=prod)))
        lt1 <- t1[, length(t)]
        res1 <- lapply(seq.int(lt1), FUN=
                       function(i) {
                           res1 <- diag(r1[i, ])
                           res1[lower.tri(res1)] <- r2[i, ]
                           res2 <- matrix(res1, ncol=ncol(res1), byrow=TRUE)
                           res1[upper.tri(res1)] <- res2[upper.tri(res2)]
                           return(res1)
                       })
        res1 <- as.array(unlist(res1))
        dim(res1) <- c(g1, g1, lt1)
        if (any(is.nan(res1[, , lt1]))) {
            res1[, , lt1][which(is.nan(res1[, , lt1]))] <- 0
        }
        dimnames(res1) <- list(x[, unique(cg)], x[, unique(cg)], t1[, t])
    }
    class(res1) <- c("COV", class(res1))
    data.table::setattr(x, "COV", res1)
    return(attr(x, "COV"))
}
#' @rdname COV
#' @method COV stratTen
#' @aliases COV.stratTen
#' @export
#' 
COV.stratTen <- function(x, ..., reCalc=FALSE){
    lapply(x, FUN=COV, reCalc=reCalc)
    lapply(x, attr, "COV")
}
#' @rdname COV
#' @method COV numeric
#' @aliases COV.numeric
#' @export
#' @examples
#' ## example of numeric method
#' ## Three covariate groups
#' ## K&M. Example 7.4, pg 212.
#' data("bmt", package="KMsurv")
#' b1 <- asWide(ten(Surv(time=t2, event=d3) ~ group, data=bmt))
#' rowSums(b1[, COV(x=e, n=n, ncg=matrix(data=c(n_1, n_2, n_3), ncol=3))], dims=2)
#' 
COV.numeric <- function(x, ..., n, ncg){
  stopifnot(all(sapply(list(x, n, ncg), is.numeric)))
  ## ensure all same length
  stopifnot(
    diff(range(sapply(list(x, n), length)))
    < .Machine$double.eps)
  ## no. of groups
  g1 <- ncol(ncg)
  if (is.null(g1)) g1 <- 1L
  ## if 2 groups only
  if (g1==1) {
    cov1 <- (ncg / n) * (1 - (ncg / n)) * ((n - x) / (n - 1)) * x
    return(cov1)
  }
  ## hold results
  a1 <- array(data=0,  dim=c(g1, g1, length(n)))
  ## diagonal elements
  for (i in seq_len(g1)) {
    a1[i, i, ] <-
      (ncg[, i] * (n - ncg[, i]) * x * (n - x)) / (n^2 * (n - 1))
  }
  ## off-diagonal elements
  for (j in seq_len(g1)) {
    for (k in 1:g1){
      if (j==k) next
      a1[j, k, ] <-
        - (ncg[, j] * ncg[, k] * x * (n - x)) / (n^2 * (n-1))
    }
  }
  if (any(is.nan(a1[, , length(n)]))) {
      a1[, , length(n)][which(is.nan(a1[, , length(n)]))] <- 0
  }
  dimnames(a1) <- list(1:g1, 1:g1, seq.int(length(n)))
  class(a1) <- c("COV", class(a1))
  return(a1)
}
