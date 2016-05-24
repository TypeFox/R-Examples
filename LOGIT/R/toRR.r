# Hilbe, J.M., Practical Guide to Logistic Regression
# coef to IRR following glm-poisson or glm.nb  Joseph M Hilbe 7Feb2015
#' @title   Display count model coefficient table as incidence rate ratios and associated statistics.
#' @description Following the glm or glm.nb commands, toRR() displays a table of incidence rate ratios and related statistics
#' including exponentiated model confidence intervals.
#' @aliases toRR
#' @usage toRR(object)
#' @importFrom stats vcov pnorm qnorm
#' @format  \describe{
#' \item{object}{
#' The only argument is the name of the fitted glm or glm.nb function model}
#'
#' \item{or}{incidence rate ratio of predictor}
#' \item{delta}{Model standard error using delta method}
#' \item{zscore}{z-statistic}
#' \item{pvalue}{probability-value based on normal distribution}
#' \item{exp.loci}{Exponentialed lower model confidence interval}
#' \item{exp.upci}{Expontiated upper model confidence interval}
#' }
#'
#'
#' @details toRR is a post-estimation function, following the use of glm() with the Poisson
#' or negative.binomial families, and following glm.nb().
#' @param object  name of the fitted glm function model
#' @return list
#' @note toRR must be loaded into memory in order to be effectve. As a function in LOGIT,
#' it is immediately available to a user.
#'@examples
#' library(MASS)
#' library(LOGIT)
#' data(medpar)
#' medpar$los<-as.numeric(medpar$los)
#' mypoi <- glm(los ~ white + hmo + factor(age80), family=poisson, data=medpar)
#' summary(mypoi)
#' toRR(mypoi)
#'
#' @seealso \code{\link{glm}}, \code{\link{glm.nb}}
#' @author Joseph M. Hilbe, Arizona State University, and
#' Jet Propulsion Laboratory, California Institute of technology
#'
#'
#'
#' @references Hilbe, Joseph M. (2015), Practical Guide to Logistic Regression, Chapman & Hall/CRC.
#'
#' Hilbe, Joseph M. (2014), Modeling Count Data, Cambridge University Press.
#' @keywords models
#' @export
#'
toRR <- function(object) {
     coef <- object$coef
       se <- sqrt(diag(vcov(object)))
       zscore <- coef / se
     rr <- exp(coef)
       delta <- rr * se
      pvalue <- 2*pnorm(abs(zscore),lower.tail = FALSE)
     loci <- coef - qnorm(.975) * se
     upci <- coef + qnorm(.975) * se
     rrtab <- data.frame(rr, delta, zscore, pvalue, exp(loci), exp(upci))
     round(rrtab, 4)
}
