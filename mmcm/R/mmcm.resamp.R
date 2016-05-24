#' The permuted modified maximum contrast method
#'
#' This function gives \eqn{P}-value for the permuted modified maximum
#' contrast method.
#'
#' \code{\link{mmcm.resamp}} performs the permuted modified maximum contrast
#' method that is detecting a true response pattern under the unequal sample size
#' situation.
#' 
#' \eqn{Y_{ij} (i=1, 2, \ldots; j=1, 2, \ldots, n_i)}{Y_ij (i = 1, 2, ...;
#' j = 1, 2, ..., n_i)} is an observed response for \eqn{j}-th individual in
#' \eqn{i}-th group.
#' 
#' \eqn{\bm{C}}{C} is coefficient matrix for permuted modified maximum contrast
#' statistics (\eqn{i \times k}{i x k} matrix, \eqn{i}: No. of groups, \eqn{k}:
#' No. of pattern).
#' \deqn{
#'   \bm{C}=(\bm{c}_1, \bm{c}_2, \ldots, \bm{c}_k)^{\rm{T}}
#' }{
#'   C = (c_1, c_2, ..., c_k)^T
#' }
#' \eqn{\bm{c}_k}{c_k} is coefficient vector of \eqn{k}-th pattern.
#' \deqn{
#'   \bm{c}_k=(c_{k1}, c_{k2}, \ldots, c_{ki})^{\rm{T}} \qquad (\textstyle \sum_i c_{ki}=0)
#' }{
#'   c_k = (c_k1, c_k2, ..., c_ki)^T (sum from i of c_ki = 0)
#' }
#' 
#' \eqn{M_{\max}}{M_max} is a permuted modified maximum contrast statistic.
#' \deqn{
#'   \bar{Y}_i=\frac{\sum_{j=1}^{n_i} Y_{ij}}{n_{i}},
#'   \bar{\bm{Y}}=(\bar{Y}_1, \bar{Y}_2, \ldots, \bar{Y}_i, \ldots, \bar{Y}_a)^{\rm{T}},
#'   M_{k}=\frac{\bm{c}^t_k \bar{\bm{Y}}}{\sqrt{\bm{c}^t_k \bm{c}_k}},
#' }{
#'   Ybar_i = (sum from j of Y_ij) / n_i,
#'   Ybar = (Ybar_1, Ybar_2, ..., Ybar_i, ..., Ybar_a)^T (a x 1 vector),
#'   M_k = c_k^t Ybar / (c_k^t c_k)^(1/2),
#' }
#' \deqn{
#'   M_{\max}=\max(M_1, M_2, \ldots, M_k).
#' }{
#'   M_max = max(M_1, M_2, ..., M_k).
#' }
#' 
#' Consider testing the overall null hypothesis \eqn{H_0: \mu_1=\mu_2=\ldots=\mu_i},
#' versus alternative hypotheses \eqn{H_1} for response petterns 
#' (\eqn{H_1: \mu_1<\mu_2<\ldots<\mu_i, \mu_1=\mu_2<\ldots<\mu_i,
#' \mu_1<\mu_2<\ldots=\mu_i}).
#' The \eqn{P}-value for the probability distribution of \eqn{M_{\max}}{M_max}
#' under the overall null hypothesis is
#' \deqn{
#'   P\mbox{-value}=\Pr(M_{\max}>m_{\max} \mid H_0)
#' }{
#'   P-value = Pr(M_max > m_max | H0)
#' }
#' \eqn{m_{\max}}{m_max} is observed value of statistics.
#' This function gives distribution of \eqn{M_{\max}}{M_max} by using the
#' permutation method, follow algorithm:
#' 
#' 1. Initialize counting variable: \eqn{COUNT = 0}{COUNT = 0}.
#' Input parameters: \eqn{NRESAMPMIN}{NRESAMPMIN} (minimum resampling count,
#' we set 1000), \eqn{NRESAMPMAX}{NRESAMPMAX} (maximum resampling count), and
#' \eqn{\epsilon}{epsilon} (absolute error tolerance).
#' 
#' 2. Calculate \eqn{m_{\max}}{m_max} that is the observed value of the test
#' statistic.
#' 
#' 3. Let \eqn{y_{ij}^{(r)}}{y_ij(r)} donate data, which are sampled without
#' replacement, and independently, form observed value \eqn{y_{ij}}{y_ij}.
#' Where, \eqn{(r)}{(r)} is suffix of the resampling number
#' \eqn{(r = 1,\, 2,\, \ldots)}{(r = 1, 2, ...)}.
#' 
#' 4. Calculate \eqn{m^{(r)}_{\max}}{m(r)_max} from \eqn{y_{ij}^{(r)}}{y_ij(r)}.
#' If \eqn{m^{(r)}_{\max} > m_{\max}}{m(r)_max > m_max}, then increment the
#' counting variable: \eqn{COUNT = COUNT + 1}{COUNT = COUNT + 1}. Calculate
#' approximate P-value \eqn{\hat{p}^{(r)}=COUNT/r}{hat-p(r) = COUNT / r}, and
#' the simulation standard error \eqn{SE(\hat{p}^{(r)})=\sqrt{\hat{p}^{(r)}(1-
#' \hat{p}^{(r)})/r}}{[hat-p(r) (1 - hat-p(r)) / r]^(1/2)}.
#' 
#' 5. Repeat 3--4, while \eqn{r > 1000}{r > 1000} and \eqn{3.5SE(\hat{p}^{(r)})
#' < \epsilon}{3.5 SE(hat-p(r)) < epsilon} (corresponding to 99\% confidence
#' level), or \eqn{NRESAMPMAX}{NRESAMPMAX} times. Output the approximate P-value
#' \eqn{\hat{p}^{(r)}}{hat-p(r)}.
#'
#' @param x a numeric vector of data values
#' @param g a integer vector giving the group for the corresponding elements of x
#' @param contrast a numeric contrast coefficient matrix for permuted modified
#' maximum contrast statistics
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less".
#' You can specify just the initial letter.
#' @param nsample specifies the number of resamples (default: 20000)
#' @param abseps specifies the absolute error tolerance (default: 0.001)
#' @param seed a single value, interpreted as an integer;
#' see \code{\link[base:Random]{set.seed()}} function. (default: NULL)
#' @return
#' \item{statistic}{the value of the test statistic with a name describing it.}
#' \item{p.value}{the p-value for the test.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{the type of test applied.}
#' \item{contrast}{a character string giving the names of the data.}
#' \item{contrast.index}{a suffix of coefficient vector of the \eqn{k}th pattern
#' that gives permuted modified maximum contrast statistics (row number of the
#' coefficient matrix).}
#' \item{error}{estimated absolute error and,}
#' \item{msg}{status messages.}
#' @references
#' Nagashima, K., Sato, Y., Hamada, C. (2011).
#' A modified maximum contrast method for unequal sample sizes in
#' pharmacogenomic studies
#' \emph{Stat Appl Genet Mol Biol.} \strong{10}(1): Article 41.
#' \url{http://dx.doi.org/10.2202/1544-6115.1560}
#' 
#' Sato, Y., Laird, N.M., Nagashima, K., et al. (2009).
#' A new statistical screening approach for finding pharmacokinetics-related
#' genes in genome-wide studies.
#' \emph{Pharmacogenomics J.} \strong{9}(2): 137--146.
#' \url{http://www.ncbi.nlm.nih.gov/pubmed/19104505}
#' @seealso
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
#'   main  = "Dominant model (sample data)",
#'   xlab  = "Genotype", ylab="PK parameter"
#' )
#' 
#' # coefficient matrix
#' # c_1: additive, c_2: recessive, c_3: dominant
#' contrast <- rbind(
#'   c(-1, 0, 1), c(-2, 1, 1), c(-1, -1, 2)
#' )
#' y <- mmcm.resamp(x, g, contrast, nsample = 20000, abseps = 0.01, seed = 5784324)
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
#' 
#' mmcmtapply <- function(r) {
#'   mmcm.resamp(
#'     xx$x[xx$pos==r[1]], xx$g[xx$pos==r[1]],
#'     contrast, nsample = 10000, abseps = 0.01, seed = 5784324+as.numeric(r[1])
#'   )
#' }
#' y <- tapply(xx$pos, xx$pos, mmcmtapply)
#' yy <- data.frame(
#'   Pos      = as.vector(names(y)),
#'   Pval     = as.vector(sapply(y, "[[", 3)),
#'   Pattern  = as.vector(sapply(y, "[[", 7)),
#'   MC_Error = as.vector(sapply(y, "[[", 9))
#' )
#' yy
#' @keywords htest
#' @importFrom stats var
#' @export
mmcm.resamp <- function(x, g, contrast, alternative = c("two.sided", "less", "greater"),
  nsample = 20000, abseps = 0.001, seed = NULL) {
  
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
  
  METHOD <- "Permuted modified maximum contrast method"
  
  p          <- length(unique(g))
  m          <- nrow(contrast)
  df         <- length(g) - p
  pooled     <- (tapply(x, g, length) - 1) * tapply(x, g, var)
  pooled     <- t(rep(1, p)) %*% pooled / df
  D          <- diag(1 / as.vector(tapply(x, g, length)))
  CDC        <- contrast %*% D %*% t(contrast)
  
  CtC        <- contrast %*% t(contrast)
  CtCMATRIX  <- matrix(rep(1, m * m), ncol = m) * diag(CtC)
  Rs         <- CDC / (sqrt(CtCMATRIX) * sqrt(t(CtCMATRIX)))
  
  STATISTICS <- switch(
    alternative,
    less      = (contrast %*% tapply(x, g, mean) / sqrt(diag(CtC))),
    greater   = (contrast %*% tapply(x, g, mean) / sqrt(diag(CtC))),
    two.sided = abs(contrast %*% tapply(x, g, mean) / sqrt(diag(CtC)))
  )
  STATISTIC  <- switch(
    alternative,
    less      = min(STATISTICS),
    greater   = max(STATISTICS),
    two.sided = max(STATISTICS)
  )
  IMAXCONT <- (1:m)[STATISTICS==STATISTIC]
  NMAXCONT <- contrast[IMAXCONT,]
  
  nalternative <- switch(
    alternative,
    less      = 1,
    greater   = 2,
    two.sided = 3
  )
  if (!is.null(seed)) {
    set.seed(seed)
  }
  RESAMP <- .C(
    "mmcm_rwrap",
    as.double(x),
    as.double(g),
    as.double(as.vector(t(contrast))),
    as.integer(nsample),
    as.integer(ncol(contrast)),
    as.integer(nrow(contrast)),
    as.integer(length(x)),
    as.double(abseps),
    as.integer(nalternative),
    pval=double(1),
    error=double(1)
  )
  
  if (RESAMP$error < abseps) {
    msg <- "Normal Completion"
  } else {
    msg <- warning("Completion with error > abseps")
  }
  
  PVAL  <- RESAMP$pval
  ERROR <- RESAMP$error
  MSG   <- msg
  
  if (length((1:m)[STATISTICS == STATISTIC]) != 1) {
    MAXCONT <- warning("More than 2 contrast coefficient vectors were selected")
  } else {
    MAXCONT <- "("
    for(i in 1:p) {
      if (i == p) {
        MAXCONT <- paste(MAXCONT, NMAXCONT[i], ")", sep = "")
      } else {
        MAXCONT <- paste(MAXCONT, NMAXCONT[i], ", ", sep = "")
      }
    }
  }
  
  names(STATISTIC) <- "Permuted modified maximum contrast statistic"
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
    msg            = MSG
    ),
    class          = "mmcm"
  )
  return(RVAL)
  
}
