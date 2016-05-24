#' Information Criteria Confidence Intervals
#'
#' Calculate confidence intervals of AIC and BIC for non-nested models.
#'
#' Functionality is currently available for models of classes
#' \code{lm}, \code{glm}, \code{glm.nb}, \code{clm}, \code{hurdle}, \code{zeroinfl}, \code{mlogit}, \code{nls}, \code{polr}, \code{rlm}, and \code{lavaan}.
#'
#' Users should take care to ensure that the two models have
#' the same dependent variable (or, for lavaan objects, identical
#' modeled variables), with observations ordered identically within
#' each model object.  Assuming the same data matrix is used to fit each model,
#' observation ordering should generally be identical.  There are currently
#' no checks for this, however.
#'
#' Note: if models are nested or if the "variance test" from
#' \code{vuongtest()} indicates models are indistinguishable, then the
#' intervals returned from \code{icci()} will be incorrect.
#'
#' @param object1 a model object
#' @param object2 a model object
#' @param conf.level confidence level of the interval
#'
#' @author Ed Merkle and Dongjun You
#'
#' @return an object of class \code{icci} containing test results.
#'
#' @references
#'
#' Vuong, Q. H. (1989).  Likelihood ratio tests for model selection and non-nested hypotheses.  \emph{Econometrica, 57}, 307-333.
#'
#' Merkle, E. C., You, D., & Preacher, K. (2014). Testing non-nested structural equation models.  \emph{Manuscript under review}.
#'
#' @examples
#' \dontrun{
#' ## Count regression comparisons
#' require(MASS)
#' house1 <- glm(Freq ~ Infl + Type + Cont, family=poisson, data=housing)
#' house2 <- glm(Freq ~ Infl + Sat, family=poisson, data=housing)
#'
#' ## CI for BIC
#' icci(house2, house1)
#'
#' ## Further comparisons to hurdle, zero-inflated models
#' require(pscl)
#' bio1 <- glm(art ~ fem + mar + phd + ment, family=poisson, data=bioChemists)
#' bio2 <- hurdle(art ~ fem + mar + phd + ment, data=bioChemists)
#' bio3 <- zeroinfl(art ~ fem + mar + phd + ment, data=bioChemists)
#' icci(bio2, bio1)
#' icci(bio3, bio1)
#' icci(bio3, bio2)
#'
#' ## Latent variable model comparisons
#' require(lavaan)
#' HS.model <- 'visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#' fit1 <- cfa(HS.model, data=HolzingerSwineford1939)
#' fit2 <- cfa(HS.model, data=HolzingerSwineford1939, group="school")
#' icci(fit1, fit2)
#' }
#'
#' @importFrom stats AIC var qnorm
#' @export
icci <- function(object1, object2, conf.level=.95) {
  classA <- class(object1)[1L]
  classB <- class(object2)[1L]
  callA <- if (isS4(object1)) object1@call else object1$call
  callB <- if (isS4(object2)) object2@call else object2$call

  llA <- llcont(object1)
  llB <- llcont(object2)

  ## Eq (4.2)
  n <- NROW(llA)
  omega.hat.2 <- (n-1)/n * var(llA - llB)

  ## BIC is computed like this because hurdle, zeroinfl, mlogit
  ## don't have an nobs() method
  bicA <- AIC(object1, k = log(NROW(estfun(object1))))
  bicB <- AIC(object2, k = log(NROW(estfun(object2))))

  aicA <- AIC(object1)
  aicB <- AIC(object2)

  bicdiff <- bicA - bicB
  aicdiff <- aicA - aicB
  alpha <- 1 - conf.level

  ## BIC CI
  BICci <- bicdiff + qnorm(c(alpha/2,(1-alpha/2)))*sqrt(n * 4 * omega.hat.2)
  ## AIC CI
  AICci <- aicdiff + qnorm(c(alpha/2,(1-alpha/2)))*sqrt(n * 4 * omega.hat.2)

  rval <- list(class = list(class1=classA, class2=classB),
               call = list(call1=callA, call2=callB),
               BIC = list(BIC1=bicA, BIC2=bicB),
               BICci = BICci,
               AIC = list(AIC1=aicA, AIC2=aicB),
               AICci = AICci,
               confLevel = conf.level)
  class(rval) <- "icci"
  return(rval)
}

################################################################
## print method for icci
################################################################
#' @method print icci
#' @export
print.icci <- function(x, ...) {
  cat("\nModel 1 \n")
  cat(" Class:", x$class$class1, "\n")
  ## a char vector with each element of length 'width.cutoff'
  model1call <- deparse(x$call$call1)
  cat(" Call: ", model1call[1], if (length(model1call) > 1) "...\n" else "\n", sep="")
  cat(" AIC:", formatC(x$AIC$AIC1, digits=3L, format="f"), "\n")
  cat(" BIC:", formatC(x$BIC$BIC1, digits=3L, format="f"), "\n")
  cat("\nModel 2 \n")
  cat(" Class:", x$class$class2, "\n")
  ## a char vector with each element of length 'width.cutoff'
  model2call <- deparse(x$call$call2)
  cat(" Call: ", model2call[1], if (length(model2call) > 1) "...\n" else "\n", sep="")
  cat(" AIC:", formatC(x$AIC$AIC2, digits=3L, format="f"), "\n")
  cat(" BIC:", formatC(x$BIC$BIC2, digits=3L, format="f"), "\n\n")

  cat(x$confLevel * 100,
      "% Confidence Interval of AIC difference (AICdiff = AIC1 - AIC2) \n", sep="")
  cat("  ", formatC(x$AICci[1], digits=3L, format="f"), " < ", "AICdiff",
      " < ", formatC(x$AICci[2], digits=3L, format="f"), "\n\n", sep="")

  cat(x$confLevel * 100,
      "% Confidence Interval of BIC difference (BICdiff = BIC1 - BIC2) \n", sep="")
  cat("  ", formatC(x$BICci[1], digits=3L, format="f"), " < ", "BICdiff",
      " < ", formatC(x$BICci[2], digits=3L, format="f"), "\n", sep="")
}
