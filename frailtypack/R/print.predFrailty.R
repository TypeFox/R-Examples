
print.predFrailty <- function(x, digits = 3, ...) 
{
# 	cl <- x$call
# 	if (!is.null(cl)){
# 		cat("\n")
# 		cat("--------- Model ---------\n")
# 		dput(cl)
# 		cat("\n")
# 	}
	if(class(x)!="predFrailty"){
		stop("The object x must be a class predFrailty.")
	}else{
		if (!is.null(cl <- x$call)){
			cat("Call:\n")
			dput(cl)
			cat("\n")
		}
		print(x$pred)
	}
}

