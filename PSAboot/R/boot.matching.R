#' Matching package implementation for bootstrapping.
#' 
#' @param estimand character string for estimand, either ATE, ATT, or ATC. See
#'        \code{\link{Match}} for more details. 
#' @param ... other parameters passed to \code{\link{Match}}.
#' @inheritParams boot.strata
#' @export
boot.matching <- function(Tr, Y, X, X.trans, formu, estimand='ATE', ...) {
	formu <- update.formula(formu, 'treat ~ .')
	ps <- fitted(glm(formu, data=cbind(treat=Tr, X), family='binomial'))
	mr <- Matching::Match(Y=Y, Tr=Tr, X=ps, estimand=estimand, ...)
	ttest <- t.test(Y[mr$index.treated], Y[mr$index.control], paired=TRUE)
	return(list(
		summary=c(estimate=unname(ttest$estimate),
				  ci.min=ttest$conf.int[1],
				  ci.max=ttest$conf.int[2],
				  t=unname(ttest$statistic),
				  p=ttest$p.value ),
		details=list(Match=mr, t.test=ttest),
		balance=balance.matching(mr$index.treated, mr$index.control, X.trans)
	))
}

