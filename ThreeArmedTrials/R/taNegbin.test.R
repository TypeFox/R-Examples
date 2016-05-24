#' @useDynLib ThreeArmedTrials
#' @import stats
#' @import MASS

#' @title Statistical test for three-armed clinical trials with negative binomial distributed endpoints.
#' @description Wald-type test for superiority/non-inferiority of the experimental treatment versus reference treatment with respect to placebo.
#' @details The hypothesis \eqn{(\lambda_P - \lambda_E)/(\lambda_P - \lambda_R) \le \Delta} is tested against the alternative \eqn{(\lambda_P - \lambda_E)/(\lambda_P - \lambda_R) > \Delta}. 
#' \eqn{\lambda_E}, \eqn{\lambda_R}, \eqn{\lambda_P} are the rates of the experimental treatment (\code{rateExp}), the reference treatment (\code{rateRef}), and the placebo group (\code{ratePla}), respectively.
#' The margin \emph{Delta}, i.e. \eqn{\Delta} in the formulas above, is between 0 and 1 for testing non-inferiority and larger than 1 for testing superiority. 
#' The parametrisation of the underlying negative binomial distributions is chosen such that 
#' a negative binomial distribution of rate \eqn{\lambda} and shape parameter \eqn{\phi} has variance \eqn{\lambda(1+\phi \lambda)}.
#' The shape parameter \eqn{\phi} is assumed to be the same among the groups. 
#' @param xExp A (non-empty) numeric vector of data values coming from the experimental treatment group.
#' @param xRef A (non-empty) numeric vector of data values coming from the reference treatment group.
#' @param xPla A (non-empty) numeric vector of data values coming from the placebo group.
#' @param Delta A numeric value specifying the non-inferiority or superiority margin. Is between 0 and 1 in case of non-inferiority and larger than 1 in case of superiority.
#' @param method A character string determing how the variance for the Wald-type test statistic is estimated, must be \emph{RML}, \emph{ML}, or \emph{SampleVariance}.
#' @return A list with class "htest" containing the following components:
#' \item{statistic}{The value of the Wald-type test statistic.}
#' \item{p.value}{The p-value for the Wald-type test.}
#' \item{method}{A character string indicating what type of Wald-type-test was performed.}
#' \item{estimate}{The estimated rates for each of the group as well as the maximum-likelihood estimator for the shape parameter.}
#' \item{sample.size}{The total number of data points used for the Wald-type test.}
#' @examples 
#' # Negative binomially distributed endpoints
#' # Test for non-inferiority test. lambda_P=8, lambda_R = 4, lambda_E = 5, and phi = 1
#' # Delta = (lambda_P-lambda_E)/(lambda_P-lambda_R)
#' xExp <-rnbinom(60, mu=5, size=1)
#' xRef <-rnbinom(40, mu=4, size=1)
#' xPla <-rnbinom(40, mu=8, size=1)
#' Delta <- (8-5)/(8-4)
#' taNegbin.test(xExp, xRef, xPla, Delta, method = 'RML')
#' taNegbin.test(xExp, xRef, xPla, Delta, method = 'ML')
#' taNegbin.test(xExp, xRef, xPla, Delta, method = 'SampleVariance')
#' 
#' # Test for superiority test. lambda_P=8, lambda_R = 5, lambda_E = 4, and phi = 1
#' # Delta = (lambda_P-lambda_E)/(lambda_P-lambda_R)
#' xExp <-rnbinom(60, mu=5, size=1)
#' xRef <-rnbinom(40, mu=4, size=1)
#' xPla <-rnbinom(40, mu=8, size=1)
#' Delta <- (8-5)/(8-4)
#' taNegbin.test(xExp, xRef, xPla, Delta, method = 'RML')
#' taNegbin.test(xExp, xRef, xPla, Delta, method = 'ML')
#' taNegbin.test(xExp, xRef, xPla, Delta, method = 'SampleVariance')
#' @references Muetze T et al. 2015. \emph{Statistical inference for three-arm trials with negative binomially distributed endpoints.} (Submitted.)
#' @seealso \code{\link{power.taNegbin.test}}
#' @export
#' @keywords test NegativeBinomial
taNegbin.test <- function(xExp, xRef, xPla, Delta, method = c('RML', 'ML', 'SampleVariance')){
  
  # To-Do: Alternative 'smaller'
  
  method <- match.arg(method)
  data.name <- paste(c(deparse(substitute(xExp)), ', ', deparse(substitute(xRef)), ', and ', deparse(substitute(xPla)), collapse=''))
                     
  xExp <- xExp[!is.na(xExp)]
  xRef <- xRef[!is.na(xRef)]
  xPla <- xPla[!is.na(xPla)]
  
  # Sample size allocation
  nExp <- length(xExp)
  nRef <- length(xRef)
  nPla <- length(xPla)
  n <- nExp + nRef + nPla
  wExp <- nExp / n
  wRef <- nRef / n
  wPla <- nPla / n
  
  # Rate estimators
  rateExp <- mean(xExp)
  rateRef <- mean(xRef)
  ratePla <- mean(xPla)
  
  # Shape parameter ML estimator
  shapeMLE <- 1/ .C("newton_Shape", zufallszahlen=as.integer(c(xExp, xRef, xPla)), mean1 = as.double(rateExp), mean2 = as.double(rateRef), mean3 = as.double(ratePla), n1 = as.integer(nExp), n2 = as.integer(nRef), n3 = as.integer(nPla), theta_out = as.double(numeric(1)) )$theta_out

  # Variance estimator
  varWaldTest <- taNegbinVarEst(xExp = xExp, xRef = xRef, xPla = xPla, Delta = Delta, method = method)
  
  # Test statistic and p-value
  TeststatMme <- sqrt(n) * (  rateExp - Delta * rateRef + (Delta-1) * ratePla ) / sqrt(varWaldTest)
  pvalueMme <- pnorm(TeststatMme)
  
  
  # Output
  estimate <- c(rateExp, rateRef, ratePla, shapeMLE)
  names(estimate) <- c('Rate Exp', 'Rate Ref', 'Rate Pla', 'Shape')
  names(TeststatMme) <- c('T')
  switch(method,
          SampleVariance = (methodText <- 'Wald-type test with sample variance based variance estimator'),
          RML= (methodText <- 'Wald-type test with restricted maximum-likelihood variance estimator'),
          ML= (methodText <- 'Wald-type test with maximum-likelihood variance estimator')
  )
  structure(list(statistic = TeststatMme,
                 p.value = pvalueMme,
                 method = methodText,
                 data.name = data.name,
                 estimate = estimate,
                 sample.size = n),
            class = "htest")
}
