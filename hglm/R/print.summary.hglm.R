`print.summary.hglm` <-
	function(x, digits = 4, print.ranef = FALSE, ...) {

x$nRand <- cumsum(x$RandC)
cat("Call: \n")
print(x$call)
cat('\n----------\n')
cat("MEAN MODEL\n")
cat('----------\n')
cat("\n")
cat("Summary of the fixed effects estimates:\n")
cat("\n")
printCoefmat(x$FixCoefMat, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
cat("Note: P-values are based on", x$devdf, "degrees of freedom\n")
if (!is.null(x$RandCoefMat)) {
	if (length(x$RandC) == 1) {
		cat("\n")
		cat("Summary of the random effects estimates:\n")
		cat("\n")
		if (nrow(x$RandCoefMat) <= 5) {
			print(round(x$RandCoefMat, digits))
		} else if (print.ranef) {
			print(round(x$RandCoefMat, digits))
		} else {
			print(round(x$RandCoefMat[1:3,], digits))
			cat('...\n')
			cat('NOTE: to show all the random effects, use print(summary(hglm.object), print.ranef = TRUE).\n')
		}
	} else {
		for (J in 1:length(x$RandC)) {
			cat("\n")
			cat("Summary of the random effects estimates:\n")
			cat("\n")
			if (nrow(x$RandCoefMat[[J]]) <= 5) {
				print(round(x$RandCoefMat[[J]], digits))
			} else if (print.ranef) {
				print(round(x$RandCoefMat[[J]], digits))
			} else {
				print(round(x$RandCoefMat[[J]][1:3,], digits))
				cat('...\n')
				cat('NOTE: to show all the random effects, use print(summary(hglm.object), print.ranef = TRUE).\n')
			}
		}
	}
}
cat('\n')
cat("----------------\n")
cat("DISPERSION MODEL\n")
cat("----------------\n")
cat("\n")
if (x$Method == "REML"){
    cat("-2*Log(Profile-Likelihood) =", round(-2*x$ProfLogLik, digits), "\n")
} else {
    cat("NOTE: h-likelihood estimates through EQL can be biased.\n")
}
if (!is.null(x$varFix)) {
	if (is.null(x$SummVC1)) {
		cat("\n") 
		cat("Dispersion parameter held Constant at:", x$varFix, '\n')
	} else {
		cat("\n") 
		cat("Dispersion parameter for the mean model:\n")
		print(x$varFix)
	}
} 
if (!is.null(x$SummVC1)) {
	cat("\n")
	cat("Model estimates for the dispersion term:\n")
	cat("\nLink =", x$link.disp, "\n")
	cat("\nEffects:\n")
	print(round(x$SummVC1, digits))
	cat("\n")
	cat("Dispersion = 1 is used in Gamma model on deviances to calculate the standard error(s).")
}
cat('\n')
if (!is.null(x$varRanef)) {
	cat("\n")
	cat("Dispersion parameter for the random effects:\n")
	print(x$varRanef, digits = digits)
	cat("\n")
	cat("Dispersion model for the random effects:\n")
	cat("\nLink = log\n")
	if (length(x$RandC) == 1) {
		cat("\nEffects:\n")
		cat(names(x$SummVC2), '\n')
		if (is.null(x$CAR.tau) & is.null(x$SAR.tau)) {
			print(round(x$SummVC2[[1]], digits)) 
		} else if (!is.null(x$CAR.tau)) {
			rownames(x$SummVC2[[1]]) <- c('1/CAR.tau', '-CAR.rho/CAR.tau')
			print(round(x$SummVC2[[1]], digits))
			cat('CAR.tau (estimated spatial variance component):', x$CAR.tau, '\n')
			cat('CAR.rho (estimated spatial correlation):', x$CAR.rho, '\n')
		} else {
			rownames(x$SummVC2[[1]]) <- c('1/sqrt(SAR.tau)', '-SAR.rho/sqrt(SAR.tau)')
			print(round(x$SummVC2[[1]], digits))
			cat('SAR.tau (estimated spatial variance component):', x$SAR.tau, '\n')
			cat('SAR.rho (estimated spatial correlation):', x$SAR.rho, '\n')
		}
		cat("\n")
	} else {
		ranefnames <- names(x$SummVC2)
		cat("\nEffects:\n")
		for (J in 1:length(x$RandC)) {
			cat(ranefnames[J], '\n')
			if (is.null(x$CAR.tau) & is.null(x$SAR.tau)) {
				print(round(x$SummVC2[[J]], digits)) 
			} else if (!is.null(x$CAR.tau)) {
				if (x$call.rand.family[[J]]$family == 'CAR') {
					rownames(x$SummVC2[[J]]) <- c('1/CAR.tau', '-CAR.rho/CAR.tau')
					print(round(x$SummVC2[[J]], digits))
					cat('CAR.tau (estimated spatial variance component):', x$CAR.tau, '\n')
					cat('CAR.rho (estimated spatial correlation):', x$CAR.rho, '\n')
				} else {
					print(round(x$SummVC2[[J]], digits)) 
				}
			} else {
				if (x$call.rand.family[[J]]$family == 'SAR') {
					rownames(x$SummVC2[[J]]) <- c('1/sqrt(SAR.tau)', '-SAR.rho/sqrt(SAR.tau)')
					print(round(x$SummVC2[[J]], digits))
					cat('SAR.tau (estimated spatial variance component):', x$SAR.tau, '\n')
					cat('SAR.rho (estimated spatial correlation):', x$SAR.rho, '\n')
				} else {
					print(round(x$SummVC2[[J]], digits)) 
				}
			}
			cat("\n")
		}
	}
	cat("Dispersion = 1 is used in Gamma model on deviances to calculate the standard error(s).\n")
}
cat("\n")
if (!is.null(x$likelihood)) {
	cat("---------------\n")
	cat("LOG-LIKELIHOODS\n")
	cat("---------------\n")
	cat("\n")
	cat("h-likelihood:", x$likelihood$hlik, "\n")
	cat("Adjusted profile likelihood", "\n") 
	cat("   Profiled over random effects:", x$likelihood$pvh, "\n")
	cat("   Profiled over fixed and random effects:", x$likelihood$pbvh, "\n")
	cat("Conditional AIC:", x$likelihood$cAIC, "\n")
	cat("\n")
}
cat(x$Method, "estimation", x$converge, "in", x$iter, "iterations.\n")

if (!is.null(x$bad)) {
	cat('\n!! Observation', x$bad, 'is too influential! Estimates are likely unreliable !!')
}

}

