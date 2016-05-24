`print.hglm` <-
	function(x, print.ranef = FALSE, ...) {

x$nRand <- cumsum(x$RandC)
cat("Call: \n")
print(x$call)
cat("\n---------------------------")
cat("\nEstimates of the mean model\n")
cat("---------------------------\n")
cat("\n")
cat("Fixed effects:\n")
print(x$fixef)
cat('\n')
if (length(x$RandC) == 1) {
	if (length(x$ranef) <= 5) {
		cat("Random effects:\n")
		print(x$ranef)
	} else if (print.ranef) {
		cat("Random effects:\n")
		print(x$ranef)
	} else {
		cat("Random effects:\n")
		print(x$ranef[1:3])
		cat('...\n')
		print(x$ranef[(x$nRand[1] - 1):x$nRand[1]])
		cat('NOTE: to show all the random effects estimates, use print(hglm.object, print.ranef = TRUE).\n')
	}
	cat('\n')
} else {
	if (length(x$ranef[1:x$nRand[1]]) <= 5) {
		cat("Random effects:\n")
		print(x$ranef[1:x$nRand[1]])
	} else if (print.ranef) {
		cat("Random effects:\n")
		print(x$ranef[1:x$nRand[1]])
	} else {
		cat("Random effects:\n")
		print(x$ranef[1:3])
		cat('...\n')
		print(x$ranef[(x$nRand[1] - 1):x$nRand[1]])
		cat('NOTE: to show all the random effects estimates, use print(hglm.object, print.ranef = TRUE).\n')
	}
	cat('\n')
	for (J in 2:length(x$RandC)) {
		if (length(x$ranef[(x$nRand[J - 1] + 1):x$nRand[J]]) <= 5) {
			cat("Random effects:\n")
			print(x$ranef[(x$nRand[J - 1] + 1):x$nRand[J]])
		} else if (print.ranef) {
			cat("Random effects:\n")
			print(x$ranef[(x$nRand[J - 1] + 1):x$nRand[J]])
		} else {
			cat("Random effects:\n")
			print(x$ranef[(x$nRand[J - 1] + 1):(x$nRand[J - 1] + 3)])
			cat('...\n')
			print(x$ranef[(x$nRand[J] - 1):x$nRand[J]])
			cat('NOTE: to show all the random effects estimates, use print(hglm.object, print.ranef = TRUE).\n')
		}
		cat('\n')
	}
}
if (!is.null(x$varFix)) {
	cat("Dispersion parameter for the mean model:", x$varFix, '\n')
} else {
	cat("------------------------------------------")
	cat("\nEstimates of the residual dispersion model\n")
	cat("------------------------------------------\n")
	cat("\nLink =", x$link.rand.disp, "\n")
	cat("\nEffects:\n")
	print(x$SummVC1[,1])
}
if (!is.null(x$varRanef)) {
	cat("\nDispersion parameter for the random effects:", x$varRanef, '\n')
} else {
	cat("\n------------------------------------------------")
	cat("\nEstimates of the random effects dispersion model\n")
	cat("------------------------------------------------\n")
	cat("\nLink =", x$link.rand.disp, "\n")
	cat("\nEffects:\n")
	for (i in 1:length(x$SummVC2)) {
		print(x$SummVC2[[i]][,1])
	}
}
cat(paste("\nEstimation", x$Converge, "in", x$iter, "iterations\n"))

if (!is.null(x$bad)) {
	cat('\n!! Observation', x$bad, 'is too influential! Estimates are likely unreliable !!')
}

}

