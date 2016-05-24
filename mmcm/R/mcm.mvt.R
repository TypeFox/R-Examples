#' The maximum contrast method by using the randomized quasi-Monte Carlo method
#'
#' This function gives \eqn{P}-value for the maximum contrast statistics by using
#' randomized quasi-Monte Carlo method from \code{\link[mvtnorm:pmvt]{pmvt}}
#' function of package \code{mvtnorm}.
#'
#' \code{\link{mcm.mvt}} performs the maximum contrast method that is detecting
#' a true response pattern.
#' 
#' \eqn{Y_{ij} (i=1, 2, \ldots; j=1, 2, \ldots, n_i)}{Y_ij (i = 1, 2, ...;
#' j = 1, 2, ..., n_i)} is an observed response for \eqn{j}-th individual in
#' \eqn{i}-th group.
#' 
#' \eqn{\bm{C}}{C} is coefficient matrix for the maximum contrast statistics
#' (\eqn{i \times k}{i x k} matrix, \eqn{i}: No. of groups, \eqn{k}: No. of pattern).
#' \deqn{
#'   \bm{C}=(\bm{c}_1, \bm{c}_2, \ldots, \bm{c}_k)^{\rm{T}}
#' }{
#'   C = (c_1, c_2, ..., c_k)^T
#' }
#' \eqn{\bm{c}_k}{c_k} is coefficient vector of \eqn{k}th pattern.
#' \deqn{
#'   \bm{c}_k=(c_{k1}, c_{k2}, \ldots, c_{ki})^{\rm{T}} \qquad (\textstyle \sum_i c_{ki}=0)
#' }{
#'   c_k = (c_k1, c_k2, ..., c_ki)^{\rm{T}} (sum from i of c_ki = 0)
#' }
#' 
#' \eqn{T_{\max}}{T_max} is the maximum contrast statistic.
#' \deqn{
#'   \bar{Y}_i=\frac{\sum_{j=1}^{n_i} Y_{ij}}{n_{i}},
#'   \bar{\bm{Y}}=(\bar{Y}_1, \bar{Y}_2, \ldots, \bar{Y}_i, \ldots, \bar{Y}_a)^{\rm{T}},
#' }{
#'   Ybar_i = (sum from j of Y_ij) / n_i,
#'   Ybar = (Ybar_1 Ybar_2 ... Ybar_i ... Ybar_a)^T (a x 1 vector),
#' }
#' \deqn{
#'   \bm{D}=diag(n_1, n_2, \ldots, n_i, \ldots, n_a),
#'   V=\frac{1}{\gamma}\sum_{j=1}^{n_i}\sum_{i=1}^{a} (Y_{ij}-\bar{Y}_i)^2,
#' }{
#'   D = diag(n_1, n_2, ..., n_i, ..., n_a) (a x a matrix),
#'   V = 1/gamma * sum_{j=1}^{n_i} sum_{i=1}^{a} (Y_ij-Ybar_i)^2,
#' }
#' \deqn{
#'   \gamma=\sum_{i=1}^{a} (n_i-1),
#'   T_{k}=\frac{\bm{c}^t_k \bar{\bm{Y}}}{\sqrt{V\bm{c}^t_k \bm{D} \bm{c}_k}},
#' }{
#'   gamma = sum_{i=1}^{a} (n_i-1),
#'   T_k = c_k^t Ybar / (V c_k^t D c_k)^(1/2),
#' }
#' \deqn{
#'   T_{\max}=\max(T_1, T_2, \ldots, T_k).
#' }{
#'   T_max = max(T_1, T_2, ..., T_k).
#' }
#' 
#' Consider testing the overall null hypothesis \eqn{H_0: \mu_1=\mu_2=\ldots=\mu_i},
#' versus alternative hypotheses \eqn{H_1} for response petterns
#' (\eqn{H_1: \mu_1<\mu_2<\ldots<\mu_i,~ \mu_1=\mu_2<\ldots<\mu_i,~
#' \mu_1<\mu_2<\ldots=\mu_i}).
#' The \eqn{P}-value for the probability distribution of \eqn{T_{\max}}{T_max}
#' under the overall null hypothesis is
#' \deqn{
#'   P\mbox{-value}=\Pr(T_{\max}>t_{\max} \mid H_0)
#' }{
#'   P-value = Pr(T_max > t_max | H0)
#' }
#' \eqn{t_{\max}}{t_max} is observed value of statistics.
#' This function gives distribution of \eqn{T_{\max}}{T_max} by using randomized
#' quasi-Monte Carlo method from package \code{mvtnorm}.
#'
#' @param x a numeric vector of data values
#' @param g a integer vector giving the group for the corresponding elements of x
#' @param contrast a numeric contrast coefficient matrix for the maximum contrast
#' statistics
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less". You can specify
#' just the initial letter.
#' @param algorithm an object of class \code{\link[mvtnorm:algorithms]{GenzBretz}}
#' defining the hyper parameters of this algorithm
#' @return
#' \item{statistic}{the value of the test statistic with a name describing it.}
#' \item{p.value}{the p-value for the test.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{the type of test applied.}
#' \item{contrast}{a character string giving the names of the data.}
#' \item{contrast.index}{a suffix of coefficient vector of the \eqn{k}th pattern
#' that gives maximum contrast statistics (row number of the coefficient matrix).}
#' \item{error}{estimated absolute error and,}
#' \item{msg}{status messages.}
#' @references
#' Yoshimura, I., Wakana, A., Hamada, C. (1997).
#' A performance comparison of maximum contrast methods to detect dose dependency.
#' \emph{Drug Information J.} \strong{31}: 423--432.
#' @seealso
#' \code{\link[mvtnorm:pmvt]{pmvt}},
#' \code{\link[mvtnorm:algorithms]{GenzBretz}},
#' \code{\link{mmcm.mvt}}
#' @examples
#' ## Example 1 ##
#' #  true response pattern: dominant model c=(1, 1, -2)
#' set.seed(136885)
#' x <- c(
#'   rnorm(130, mean =  1 / 6, sd = 1),
#'   rnorm( 90, mean =  1 / 6, sd = 1),
#'   rnorm( 10, mean = -2 / 6, sd = 1)
#' )
#' g <- rep(1:3, c(130, 90, 10))
#' boxplot(
#'   x ~ g,
#'   width = c(length(g[g==1]), length(g[g==2]), length(g[g==3])),
#'   main = "Dominant model (sample data)",
#'   xlab = "Genotype",
#'   ylab = "PK parameter"
#' )
#' 
#' # coefficient matrix
#' # c_1: additive, c_2: recessive, c_3: dominant
#' contrast <- rbind(
#'   c(-1, 0, 1), c(-2, 1, 1), c(-1, -1, 2)
#' )
#' y <- mcm.mvt(x, g, contrast)
#' y
#' 
#' ## Example 2 ##
#' #  for dataframe
#' #  true response pattern: pos = 1 dominant  model c=( 1,  1, -2)
#' #                               2 additive  model c=(-1,  0,  1)
#' #                               3 recessive model c=( 2, -1, -1)
#' set.seed(3872435)
#' x <- c(
#'   rnorm(130, mean =  1 / 6, sd = 1),
#'   rnorm( 90, mean =  1 / 6, sd = 1),
#'   rnorm( 10, mean = -2 / 6, sd = 1),
#'   rnorm(130, mean = -1 / 4, sd = 1),
#'   rnorm( 90, mean =  0 / 4, sd = 1),
#'   rnorm( 10, mean =  1 / 4, sd = 1),
#'   rnorm(130, mean =  2 / 6, sd = 1),
#'   rnorm( 90, mean = -1 / 6, sd = 1),
#'   rnorm( 10, mean = -1 / 6, sd = 1)
#' )
#' g   <- rep(rep(1:3, c(130, 90, 10)), 3)
#' pos <- rep(c("rsXXXX", "rsYYYY", "rsZZZZ"), each=230)
#' xx  <- data.frame(pos = pos, x = x, g = g)
#' 
#' # coefficient matrix
#' # c_1: additive, c_2: recessive, c_3: dominant
#' contrast <- rbind(
#'   c(-1, 0, 1), c(-2, 1, 1), c(-1, -1, 2)
#' )
#' mmcmtapply <- function(r) {
#'   mcm.mvt(
#'     xx$x[xx$pos==r[1]], xx$g[xx$pos==r[1]],
#'     contrast
#'   )
#' }
#' y <- tapply(xx$pos, xx$pos, mmcmtapply)
#' yy <- data.frame(
#'   Pos       = as.vector(names(y)),
#'   Pval      = as.vector(sapply(y, "[[", 3)),
#'   Pattern   = as.vector(sapply(y, "[[", 7)),
#'   QMC_Error = as.vector(sapply(y, "[[", 9))
#' )
#' # miss-detection!
#' yy
#' @keywords htest
#' @importFrom stats var
#' @export
mcm.mvt <- function(x, g, contrast, alternative = c("two.sided", "less", "greater"), algorithm = GenzBretz()) {
  
  ####################
  # executable check
  ####################
  
  alternative <- match.arg(alternative)
  
  DNAMEX <- deparse(substitute(x))
  DNAMEG <- deparse(substitute(g))
  DNAMEC <- deparse(substitute(contrast))
  DNAME  <- paste("'", DNAMEX, "' by group '", DNAMEG,
                  "' with contrast coefficient matrix '",
                  DNAMEC, "'", sep="")
  
  if (!is.numeric(x)) {
    stop(paste(DNAMEX, "must be numeric"))
  }
  if (!is.numeric(g)) {
    stop(paste(DNAMEG, "must be numeric"))
  }
  if (!is.matrix(contrast)) {
    stop(paste(DNAMEC, "must be a matrix"))
  }
  
  x <- x[is.finite(x)]
  g <- g[is.finite(g)]
  
  if (length(x) < 1L) {
    stop(paste("not enough (finite) ", DNAMEX, "observations"))
  }
  
  if (length(x) != length(g)) {
    stop(paste(DNAME, "and", DNAMEG, "must have the same length"))
  }
  
  if (length(unique(g)) != ncol(contrast)) {
    stop(paste("nrow(", DNAMEC, ") and length(unique(", DNAMEG,
               ")) must have the same length", sep=""))
  }
  
  if (length((1:nrow(contrast))[apply(contrast, 1, sum) != rep(0, nrow(contrast))]) != 0) {
    stop("sum of contrast vector element must be 0\n")
  }
  
  ####################
  # execute mmcm
  ####################
  
  METHOD <- "Maximum contrast method"
  
  p          <- length(unique(g))
  m          <- nrow(contrast)
  df         <- length(g) - p
  pooled     <- (tapply(x, g, length) - 1) * tapply(x, g, var)
  pooled     <- t(rep(1, p)) %*% pooled / df
  D          <- diag(1 / as.vector(tapply(x, g, length)))
  CDC        <- contrast %*% D %*% t(contrast)
  if (m == 1) {
    SQRTCDCINV <- 1/sqrt(CDC)
    Rt         <- SQRTCDCINV^2 * CDC
  } else {
    SQRTCDCINV <- diag(1/sqrt(diag(CDC)))
    Rt         <- SQRTCDCINV %*% CDC %*% SQRTCDCINV
  }
  
  STATISTICS <- switch(
    alternative,
    less      = (contrast %*% tapply(x, g, mean) / sqrt(diag(CDC))) /
                (rep(1, m) %*% sqrt(pooled)),
    greater   = (contrast %*% tapply(x, g, mean) / sqrt(diag(CDC))) /
                (rep(1, m) %*% sqrt(pooled)),
    two.sided = abs(contrast %*% tapply(x, g, mean) / sqrt(diag(CDC))) /
                (rep(1, m) %*% sqrt(pooled))
  )
  STATISTIC  <- switch(
    alternative,
    less      = min(STATISTICS),
    greater   = max(STATISTICS),
    two.sided = max(STATISTICS)
  )
  IMAXCONT <- (1:m)[STATISTICS==STATISTIC]
  NMAXCONT <- contrast[IMAXCONT,]
  PVAL <- switch(
    alternative,
    less      = 1 - pmvt(lower = rep(STATISTIC, m), upper = rep(Inf, m),
                df = df, sigma = Rt, algorithm = algorithm),
    greater   = 1 - pmvt(lower = rep(-Inf, m), upper = rep(STATISTIC, m),
                df = df, sigma = Rt, algorithm = algorithm),
    two.sided = 1 - pmvt(lower = rep(-STATISTIC, m), upper = rep(STATISTIC, m),
                df = df, sigma = Rt, algorithm = algorithm)
  )
  ERROR <- attr(PVAL, "error")
  MSG   <- attr(PVAL, "msg")
  
  if (length((1:m)[STATISTICS==STATISTIC]) != 1) {
    MAXCONT <- warning("More than 2 contrast coefficient vectors were selected")
  } else {
    MAXCONT <- "("
    for(i in 1:p) {
      if (i==p) {
        MAXCONT <- paste(MAXCONT, NMAXCONT[i], ")", sep="")
      } else {
        MAXCONT <- paste(MAXCONT, NMAXCONT[i], ", ", sep="")
      }
    }
  }
  
  names(STATISTIC) <- "Maximum contrast statistic"
  names(IMAXCONT)  <- "index"
  names(MAXCONT)   <- "Maximum contrast coefficient vector"
  names(ERROR)     <- "Estimated absolute error of P-value"
  names(MSG)       <- "Status messages of P-value calculation"
  
  RVAL <- structure(list(
    statistic      = STATISTIC,
    parameter      = NULL,
    p.value        = as.numeric(PVAL), 
    alternative    = alternative,
    method         = METHOD, 
    data.name      = DNAME,
    contrast       = MAXCONT,
    contrast.index = IMAXCONT,
    error          = ERROR,
    msg            = MSG),
    class          = "mmcm"
  )
  return(RVAL)
  
}
