`print.mlds` <-
function(x, 
#x, obj of class mlds
		digits = max(3, getOption("digits") - 4), ...) {
	names(x$pscale) <- x$stimulus
	cat("\nPerceptual Scale:\n")
	print.default(format(x$pscale, digits = digits), 
		quote = FALSE, ...)
	cat("\nsigma:\n")
	print.default(format(x$sigma, digits = digits, ...), 
		quote = FALSE)
	invisible(x)
		}

`print.mlbs` <- function(x, digits = max(3, getOption("digits") - 4), ...){
	names(x$pscale) <- x$stimulus
    cat("\nPerceptual Scale:\n")
    print.default(format(x$pscale, digits = digits), quote = FALSE, 
        ...)
    cat("\nsigma:\n")
    print.default(format(x$sigma, digits = digits, ...), quote = FALSE)
    invisible(x)	
}
