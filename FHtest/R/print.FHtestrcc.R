print.FHtestrcc <-
function (x, digits = max(options()$digits - 4, 3), ...)
{
	saveopt <- options(digits = digits)
	on.exit(options(saveopt))
	if (!inherits(x, "FHtestrcc"))
                stop("Object is not of class FHtestrcc")
	cat("\n")
	writeLines(x$information)
	cat("\n")
	cat(x$data.name, sep = "\n")
	cat("\n")
	otmp <- x$obs
	etmp <- x$exp
	if (substr(x$information,2,3)=="Tr"){
	temp <- cbind(x$n, otmp, etmp, otmp - etmp)
	dimnames(temp) <- list(names(x$n), c("N", "Observed", "Expected", "O-E"))
	print(temp)
	}
	if (substr(x$information,2,3)!="Tr"){
	temp <- cbind(x$n, otmp, etmp, otmp - etmp, ((otmp - etmp)^2)/etmp,  ((otmp - etmp)^2)/diag(as.matrix(x$var)))
	dimnames(temp) <- list(names(x$n), c("N", "Observed", "Expected", "O-E", "(O-E)^2/E", "(O-E)^2/V"))
	print(temp)
	}
	if (substr(x$information,2,2)=="K"){
        cat("\nChisq= ", format(round(x$statistic, 1)), " on ",length(x$n)-1 ," degrees of freedom, p-value= ", format(signif (x$pvalue, digits)), "\n",sep="")
        }
	else
        cat("\nStatistic Z= ", format(round(x$statistic, 1)), ", p-value= ", format(signif(x$pvalue, digits)), "\n",sep="")
	cat(x$alt.phrase, sep = "\n")
	cat("\n")
        invisible(x)
}
