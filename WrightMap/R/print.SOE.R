print.SOE <-
function(x, ...) {
	scipen <- getOption("scipen")
	options(scipen = 999)
	cat("\nSummary of estimation\n\n")
	cat("Estimation method:", x$method, "\n")
	cat("Assumed population distribution:", x$distribution, "\n")
	cat("Constraint:", x$constraint, "\n\n")

	cat("Termination criteria:\n     ", x$max.iterations, "iterations\n     ", x$parameter.change, "change in parameters\n     ", 
		x$deviance.change, "change in deviance\n     ", x$max.iterations.no.improvement, "iterations without a deviance improvement\n     ", 
		x$max.steps, "Newton steps in M-step\n")

	m <- regexec("Iterations terminated because ", x$termination.reason)
	reason <- paste(unlist(regmatches(x$termination.reason, m, invert = TRUE)), collapse = "")
	cat("Estimation terminated after ", x$iterations, " iterations because ", reason, ".\n\n", sep = "")

	cat("Random number generation seed:", x$seed, "\n")
	cat(x$PV.nodes, "nodes used for drawing", x$n.plausible.values, "plausible values", "\n")
	cat(x$fit.nodes, "nodes used when computing fit", "\n")
	cat("Value for obtaining finite MLEs for zero/perfects:", x$zero.perfect.value, "\n")
	cat("\n")
	options(scipen = scipen)

}
