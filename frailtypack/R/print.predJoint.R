
print.predJoint <- function(x, digits = 3, ...) 
{
# 	cl <- x$call
# 	if (!is.null(cl)){
# 		cat("\n")
# 		cat("--------- Model ---------\n")
# 		dput(cl)
# 		cat("\n")
# 	}
	if(class(x)!="predJoint"){
		stop("The object x must be a class predJoint.")
	}else{
		if (!is.null(cl <- x$call)){
			cat("Call:\n")
			dput(cl)
			cat("\n")
		}
		if (!x$intcens) {
			cat("\n")
			cat("--------- Prediction 1 (exactly j recurrences) ---------\n")
			cat("\n")
			print(x$pred1,row.names=F,digits=digits)
			
			cat("\n")
			cat("--------- Prediction 2 (at least j recurrences) ---------\n")
			cat("\n")
			print(x$pred2,row.names=F,digits=digits)
			
			cat("\n")
			cat("--------- Prediction 3 (only parameters) ---------\n")
			cat("\n")
			print(x$pred3,row.names=F,digits=digits)
		}else{
			cat("\n")
			cat("--------- Prediction  ---------\n")
			cat("\n")
			print(x$pred2,row.names=F,digits=digits)
		}
	}
}

