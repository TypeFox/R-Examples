`print.summary.mlds` <-
function(x, 
	digits = max(3, getOption("digits") - 4), ...) {
#x, object of class summary.mlds
	cat("\nMethod:\t")
	cat(x$method)
	cat("\t\tLink:\t")
	cat(x$link)
	if (x$method == "formula") {
		cat("\n\nformula:\t")
		cat(deparse(x$formula))
		cat("\n\nparameters:\t")
		cat(format(x$par, digits = digits))
		}
	cat("\n\nPerceptual Scale:\n")
	print.default(format(x$pscale, digits = digits), 
		quote = FALSE, ...)
	cat("\nsigma:\t")
	cat(format(x$sigma, digits = digits))
	cat("\nlogLik:\t")
	cat(format(x$logLik, digits = digits)) 
	cat("\n")
	}

`print.summary.mlbs` <-
function(x, 
	digits = max(3, getOption("digits") - 4), ...) {
#x, object of class summary.mlds
	cat("\nMethod:\t")
	cat(x$method)
	cat("\t\tLink:\t")
	cat(x$link)
	if (x$method == "formula") {
		cat("\n\nformula:\t")
		cat(deparse(x$formula))
		cat("\n\nparameters:\t")
		cat(format(x$par, digits = digits))
		}
	cat("\n\nPerceptual Scale:\n")
	print.default(format(x$pscale, digits = digits), 
		quote = FALSE, ...)
	cat("\nsigma:\t")
	cat(format(x$sigma, digits = digits))
	cat("\nlogLik:\t")
	cat(format(x$logLik, digits = digits)) 
	cat("\n")
	}


