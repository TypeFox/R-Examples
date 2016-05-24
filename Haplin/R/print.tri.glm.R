"print.tri.glm"<-
function(x, ...)
{
# PRINTS THE RESULT OF A TRIAD ESTIMATION USING glm
#
##	cat("Object of class tri.glm\n")
	cat("\nNumber of triads:", round(x$ntri), "\n")
	cat("\nNumber of haplotypes:", round(x$nall), "\n")
	cat("\nParameter estimates:\n")
	print(x$result$coefficients, ...)
	invisible(x)
}
