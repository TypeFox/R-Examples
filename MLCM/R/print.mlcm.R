`print.mlcm` <-
function(x, digits = max(3, getOption("digits") - 4), ...){
	cat("\nMaximum Likelihood Conjoint Measurement\n")
	cat("\nModel:\t", switch(x$model, add = "Additive",
		ind = "Independence", full = "Full"))
	cat("\nPerceptual Scale:\n")
    print.default(format(x$pscale, digits = digits), quote = FALSE, 
        ...)
	cat("\nsigma:\n")
    print.default(format(x$sigma, digits = digits, ...), quote = FALSE)
    invisible(x)
}

