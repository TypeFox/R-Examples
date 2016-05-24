print.HWE.test <- function(x, digits){
##
## PRINTS THE RESULT OF A HWE TEST
#
if(x$fail){
	print(x$warnings)
	return(invisible())
}
.tab <- x$table
.geno <- paste(.tab$a1, .tab$a2, sep = "|")
.tab.ut <- t(round(as.matrix(.tab[,c("freq", "exp")]),1))
dimnames(.tab.ut) <- list(c("Observed", "Expected"), .geno)


if(F){
	if(!missing(digits)) {
		.dig <- .dig.overall <- digits
	}
##	else .dig.overall <- options()$digits
	else {.dig.overall <- 4
		.dig <- 2
	}
}

## if(missing(digits)) digits <- 4
digits <- 4

print(.tab.ut)
cat("Chisq: ", format(x$chisq, digits = digits), "   ", "Df: ", x$df, "   ", "P-value: ", format(x$p.value, digits = digits), "\n")

return(invisible())
}

