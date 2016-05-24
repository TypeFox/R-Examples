#' Stratification implementation for bootstrapping.
#' 
#' @param Tr vector indicating treatment assignment.
#' @param Y vector of outcome.
#' @param X matrix or data frame of covariates.
#' @param X.trans a data frame of \code{X} with factors recoded. See \code{\link{cv.trans.psa}}
#' @param formu the formula to use to estimate propensity scores. Note that the
#'        dependent varaible (i.e. treatment variable) name will be updated using
#'        the \code{Tr} vector.
#' @param nstrata number of strata to divide the propensity scores.
#' @param ... other parameters passed from \code{\link{PSAboot}}
#' @return a list with three elements:
#'         \describe{
#'         \item{\code{summary}}{a named numeric vector (with at minimum \code{estimate}, 
#'         \code{ci.min}, and \code{ci.max} but other values allowed)}
#'         \item{\code{balance}}{a named numeric vector with one element per 
#'         covariate listed in \code{X.trans} representing a balance statistic 
#'         (usually standardized effect size after adjustment)}
#'         \item{\code{details}}{an arbitrary object that contains the full results of the
#'         analysis}
#'         }
#' @export
boot.strata <- function(Tr, Y, X, X.trans, formu, nstrata=5, ...) {
	formu <- update.formula(formu, 'treat ~ .')
	ps <- fitted(glm(formu, data=cbind(treat=Tr, X), family='binomial'))
	strata <- cut(ps, quantile(ps, seq(0, 1, 1/nstrata)), include.lowest=TRUE, 
				  labels=letters[1:nstrata])
	strata.results <- psa.strata(Y=Y, Tr=Tr, strata=strata, ...)
	return(list(
		summary=c(estimate=strata.results$ATE,
				  ci.min=strata.results$CI.95[1],
				  ci.max=strata.results$CI.95[2],
				  se.wtd=strata.results$se.wtd,
				  approx.t=strata.results$approx.t),
		details=strata.results,
		balance=TriMatch::covariateBalance(X.trans, Tr, ps, strata)$effect.sizes[,'stES_adj'] ))
}
