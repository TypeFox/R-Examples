
print.Diffepoce <- function(x, digits = 3, ...) 
{
	if(class(x)!="Diffepoce"){
		stop("The object x must be a class Diffepoce.")
	}else{
# 		if (!is.null(cl <- x$call)){
# 			cat("Call:\n")
# 			dput(cl)
# 			cat("\n")
# 		}
		out <- cbind(x$pred.times,x$DEPOCE,x$TIinf,x$TIsup)
		rownames(out) <- rep(" ", length(out[,1]))
		if (x$new.data) colnames(out) <- c("pred.times","Diff MPOL","95%TIinf","95%TIsup")
		else colnames(out) <- c("pred.times","Diff CVPOL","95%TIinf","95%TIsup")
		print(out)
	}
}

