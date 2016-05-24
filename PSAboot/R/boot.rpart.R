#' Stratification using classification trees for bootstrapping.
#' 
#' @inheritParams boot.strata
#' @param minStrata minimum number of treatment or control unitis within a strata 
#'        to include that strata.
#' @export
boot.rpart <- function(Tr, Y, X, X.trans, formu, minStrata=5, ...) {
	formu <- update.formula(formu, 'treat ~ .')
	tree <- rpart::rpart(formu, data=cbind(treat=Tr, X))
	strata <- tree$where
	sizes <- reshape2::melt(table(strata, Tr))
	smallStrata <- sizes[sizes$value < minStrata,]$strata
	if(length(smallStrata) > 0) {
		rows <- !strata %in% smallStrata
		Tr <- Tr[rows]
		Y <- Y[rows]
		X <- X[rows,]
		X.trans <- X.trans[rows,]
		strata <- strata[rows]
	}
	if(length(unique(strata)) < 2) {
		stop('Classification tree (rpart) with no splits occurred.')	
	}
	strata.results <- psa.strata(Y=Y, Tr=Tr, strata=strata, ...)
	return(list(
		summary=c(estimate=strata.results$ATE,
				  ci.min=strata.results$CI.95[1],
				  ci.max=strata.results$CI.95[2],
				  se.wtd=strata.results$se.wtd,
				  approx.t=strata.results$approx.t),
		details=strata.results,
		balance=TriMatch::covariateBalance(X.trans, Tr, predict(tree), 
										   strata)$effect.sizes[,'stES_adj'] ))
}
