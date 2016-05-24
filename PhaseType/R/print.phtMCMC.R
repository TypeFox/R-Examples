print.phtMCMC <- function(x, ...) {
	cat("MCMC chain of length", x$iterations, "iterations for Phase-type inference on matrix structure:\n")
	print(x$T, quote=FALSE)
	cat("\nPriors:\n")
	for(i in 1:length(x$vars)) {
		cat(x$vars[i], " ~ Gamma(shape=", x$nu[[i]], ", scale=1/", x$zeta[[i]], ")\n", sep="")
	}
	cat("\nOriginal data consist", length(x$data), "observations\n")
	cat("\nPosterior sample mean and standard deviation for each variable, plus standard error of the mean:\n")
	print(summary(x$samples)[[1]], ...)
	cat("\nPosterior quantiles for each variable:\n")
	print(summary(x$samples)[[2]], ...)
}
