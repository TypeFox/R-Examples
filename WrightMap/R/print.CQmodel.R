print.CQmodel <-
function(x, ...) {

	catnull <- function(x, y, ...) {
		if (!is.null(x)) 
			cat(y, ...)
	}

	cat("\n")
	cat("ConQuest Output Summary:\n")
	cat("========================\n")
	catnull(x$title, x$title, "\n\n")
	catnull(x$equation, "The item model:", x$equation, "\n")
	catnull(x$nDim, x$nDim, ifelse(x$nDim == 1, "dimension", "dimensions"), "\n")
	catnull(x$participants, x$participants, "participants\n")
	catnull(x$deviance, "Deviance: ", x$deviance, " (", x$parameters, " parameters)\n\n", sep = "")
	cat("Additional information available:\n")
	catnull(x$SOE, "Summary of estimation: $SOE\n")
	catnull(x$RMP, "Response model parameter estimates: $RMP\n")
	catnull(x$reg.coef, "Regression coefficients: $reg.coef\n")
	catnull(x$variances, "Variances: $variances\n")
	catnull(x$cor.matrix, "Correlation matrix: $cor.matrix\n")
	catnull(x$cov.matrix, "Covariance matrix: $cov.matrix\n")
	catnull(x$rel.coef, "Reliabilities: $rel.coef\n")
	catnull(x$GIN, "GIN tables (thresholds): $GIN\n")
	catnull(x$GIN.deltas, "GIN tables (deltas): $GIN.deltas\n")
	catnull(x$p.est, paste(x$p.est.type, " table: $p.est\n", sep = ""))
	catnull(x$run.details, "Additional details: $run.details\n")
	cat("\n")
}
