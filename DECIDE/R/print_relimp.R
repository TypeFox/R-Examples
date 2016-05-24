`print_relimp` <-
function(dataset) {
	ri1 <- relative.importance(dataset)
	print(c("Total sample size: ", ri1$sample_size), quote=FALSE)
	print(c("Number of classes: ", ri1$no_classes), quote=FALSE)
	print(c("Size of classes: ", ri1$class_size), quote=FALSE)
	print("   ", quote=FALSE)
	print(c("Overall percentage that made the transition: ", format(ri1$percentage_overall, digits=5)), quote=FALSE)
	print(c("Percentage that made the transition for each class: ", format(ri1$percentage_class, digits=5)), quote=FALSE)
	print("   ", quote=FALSE)
	print(c("50% point for each class: ", format(ri1$fifty_point, digits=4)), quote=FALSE)
	print("   ", quote=FALSE)
	print("Parameters of logistic regression and normal distribution for each class: ", quote=FALSE)
	print(ri1$parameters, quote=FALSE)
	print("   ", quote=FALSE)

	print("Probabilities of transition: ", quote=FALSE)
	print(ri1$transition_prob, quote=FALSE)
	print("   ", quote=FALSE)

	print("Log odds of transition: ", quote=FALSE)
	print(ri1$log_odds, quote=FALSE)
	print("   ", quote=FALSE)

	print("Standard errors of log odds estimates: ", quote=FALSE)
	print(ri1$se_logodds, quote=FALSE)
	print("   ", quote=FALSE)

	print("95% confidence intervals for Log odds: ", quote=FALSE)
	print(ri1$ci_logodds, quote=FALSE)
	print("   ", quote=FALSE)

	print("Log odds ratios: ", quote=FALSE)
	print(ri1$log_oddsratios, quote=FALSE)
	print("   ", quote=FALSE)

	print("Standard errors for log odds ratios: ", quote=FALSE)
	print(ri1$se_logoddsratios, quote=FALSE)
	print("   ", quote=FALSE)

	print("95% confidence intervals for log odds ratios: ", quote=FALSE)
	print(ri1$ci_logoddsratios, quote=FALSE)
	print("   ", quote=FALSE)

	print("Odds ratios: ", quote=FALSE)
	print(ri1$oddsratios, quote=FALSE)
	print("   ", quote=FALSE)

	print("Relative importance of primary effects - 1: ", quote=FALSE)
	print(ri1$rel_imp_prim1, quote=FALSE)
	print("   ", quote=FALSE)

	print("Relative importance of primary effects - 2: ", quote=FALSE)
	print(ri1$rel_imp_prim2, quote=FALSE)
	print("   ", quote=FALSE)

	print("Relative importance of primary effects - average: ", quote=FALSE)
	print(ri1$rel_imp_prim_avg, quote=FALSE)
	print("   ", quote=FALSE)

	print("Relative importance of secondary effects - 1: ", quote=FALSE)
	print(ri1$rel_imp_sec1, quote=FALSE)
	print("   ", quote=FALSE)

	print("Relative importance of secondary effects - 2: ", quote=FALSE)
	print(ri1$rel_imp_sec2, quote=FALSE)
	print("   ", quote=FALSE)

	print("Relative importance of secondary effects - average: ", quote=FALSE)
	print(ri1$rel_imp_sec_avg, quote=FALSE)
	print("   ", quote=FALSE)

	print("Standard errors for relative importance - 1: ", quote=FALSE)
	print(ri1$se.ri.1, quote=FALSE)
	print("   ", quote=FALSE)

	print("95% confidence intervals for relative importance - 1: ", quote=FALSE)
	print(ri1$ci.ri.1, quote=FALSE)
	print("   ", quote=FALSE)

	print("Standard errors for relative importance - 2: ", quote=FALSE)
	print(ri1$se.ri.2, quote=FALSE)
	print("   ", quote=FALSE)

	print("95% confidence intervals for relative importance - 2: ", quote=FALSE)
	print(ri1$ci.ri.2, quote=FALSE)
	print("   ", quote=FALSE)

	print("Standard errors for relative importance - average: ", quote=FALSE)
	print(ri1$se.ri.avg, quote=FALSE)
	print("   ", quote=FALSE)

	print("95% confidence intervals for relative importance - average: ", quote=FALSE)
	print(ri1$ci.ri.avg, quote=FALSE)

}

