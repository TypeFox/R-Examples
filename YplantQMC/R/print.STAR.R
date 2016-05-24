#'@method print STAR
#'@S3method print STAR
print.STAR <- function(x, ...){
	cat("\nAverage displayed and projected leaf area.\n")
	cat(paste(c(rep("-",30),"\n"),collapse=""))
	cat("                         Total leaf area (m2) = ", 10^-6 * x$LA, "\n")
	if(x$sphericalSTAR){
		cat("Spherically averaged displayed leaf area (m2) = ", 10^-6 * x$DAbar, "\n")
		cat("                Spherically averaged STAR (-) = ", x$STARbar, "\n\n")
		cat("          Spherically averaged projection (-) = ", x$PALAbar, "(Should be 0.5).\n")
	} else {
		cat("             Average displayed leaf area (m2) = ", 10^-6 * x$DAbar, "\n")
		cat("                             Average STAR (-) = ", x$STARbar, "\n\n")
	}
}
