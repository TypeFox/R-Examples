# Hilbe, J.M., Practical Guide to Logistic Regression
# From coef table to OR after glm JM Hilbe 7Feb2015
#' @title  Display logistic coefficient table as odds ratios and associated statistics.
#' @description Following the glm command, toOR() displays a table of odds ratios and related statistics
#' including exponentiated model confidence intervals.
#' @aliases toOR
#' @usage toOR(object)
#' @importFrom stats vcov pnorm qnorm
#' @format  \describe{
#' \item{object}{
#' The only argument is the name of the fitted glm function model}
#' value
#'
#'  \item{or}{odds ratio of predictor}
#'  \item{delta}{Model standard error using delta method}
#'  \item{zscore}{z-statistic}
#'  \item{pvalue}{probability-value based on normal distribution}
#'  \item{exp.loci}{Exponentialed lower model confidence interval}
#'  \item{exp.upci}{Expontiated upper model confidence interval}
#' }
#'
#' @details toOR is a post-estimation function, following the use of glm().
#' @param object  name of the fitted glm function model
#' @return list
#' @note toOR must be loaded into memory in order to be effectve. As a function in LOGIT,
#' it is immediately available to a user.
#'@examples
#'  library(MASS)
#'  library(LOGIT)
#'  data(medpar)
#'  mylogit <- glm(died ~ los + white + hmo, family=binomial, data=medpar)
#'  summary(mylogit)
#'  toOR(mylogit)
#'
#' @seealso \code{\link{glm}}
#' @author Joseph M. Hilbe, Arizona State University, and Jet Propulsion Laboratory, California Institute of technology
#'
#'
#' @references Hilbe, Joseph M. (2015), Practical Guide to Logistic Regression, Chapman & Hall/CRC.
#' @keywords models
#' @export
#'
toOR <- function(object) {
     coef <- object$coef
     se <- sqrt(diag(vcov(object)))
     zscore <- coef / se
     or <- exp(coef)
     delta <- or * se
     pvalue <- 2*pnorm(abs(zscore),lower.tail=FALSE)
     loci <- coef - qnorm(.975) * se
     upci <- coef + qnorm(.975) * se
     ortab <- data.frame(or, delta, zscore, pvalue, exp(loci), exp(upci))
     round(ortab, 4)
}
