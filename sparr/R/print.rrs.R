print.rrs <- function(x, ...){
	
	if(!x$log) cat("Relative risk function.\n\n")
	else cat("Log-Relative risk function.\n\n")

	cat("--Numerator (case) density--\n")
	print.bivden(x$f)
	cat("\n--Denominator (control) density--\n")
	print.bivden(x$g)

}