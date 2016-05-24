get.chisq.stats <- function(groups, survobj) {

	# verify we have appropriate information
	if (nrow(survobj) != length(groups)) { return( rep(NA,3) ); }

	# survdiff between the groups
	survdiff.obj <- survdiff(survobj ~ groups);

	# compute the logrank P
	logrank.p <- 1 - pchisq( survdiff.obj$chisq, df = (length(survdiff.obj$n) - 1) );

	# return results
	return(
		c(
			"Chisq" = survdiff.obj$chisq,
			"DOF" = (length(survdiff.obj$n) - 1),
			"LogRank.P" = logrank.p
			)
		);
	}