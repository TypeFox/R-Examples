print.summary.fptl <-
function (x, ...) 
{
      if (!is.summary.fptl(x)) 
      	stop(paste(sQuote("x"), "is not of class", shQuote("summary.fptl")))

	Args <- formals(summary.fptl)
	labelCols = c("t[i]", "t[i]*", "tmax[i]^-", "tmax[i]", "tmax[i]^+")
	label <- format(c("Maximum slope required between two points to consider that a growing function is constant", 
					"Ratio of the global increase of the FPTL function in the growth subinterval [t[i], tmax[i]]", 
					"that should be reached to consider that it begins to increase significantly",
					"Maximum allowable distance between tmax[i] and tmax[i]^+"))
	
	cutoff <- options()$deparse.cutoff
	options(deparse.cutoff = 175)
	
	Call <- attr(x, "Call")
	FPTLCall <- attr(x, "FPTLCall")
	X0 <- paste("x0 =", format(sapply(FPTLCall, "[[", "x0"), ...))

	A <- Call[[1]][-1]
	Args[names(A)] <- A
	logic <- is.name(A$object)
	
	for (i in 1:length(x)){
		cat("\n")
		cat(X0[i])
		cat("\n\n    Interesting time instants:\n\n")
		y <- x[[i]]$instants			
		dimnames(y) <- list(rep("    ", nrow(y)), labelCols)			
		print(y, ...)

		cat("\n    Call:")
		cat("\n    ")
 		print(Call[[i]])

		if (logic){
			cat("\n    attr(", A$object, ", ", shQuote("Call"), "):", sep="")
			cat("\n    ")   		
			print(FPTLCall[[i]])		
		}
	}

	label <- paste(label, c(as.character(Args$zeroSlope), paste("10^(-", Args$p0.tol, ")", sep = ""), "", paste(Args$k, "*(tmax[i] - t[i]*)(1-FPTL(tmax[i]))", 
				sep = "")), sep = "   ")

	cat("\n\n")
	cat(label, sep= "\n")

	dp <- FPTLCall[[1]][["dp"]]
	if (is.name(dp)) cat("\n", dp, ":", sep="") else cat("\ndp:")		
	print(attr(x, "dp"))

	v <- attr(x, "vars")
    	if (!is.null(v)){
		if (logic) cat("Values of names which occur in attr(", A$object, ", ", shQuote("Call"), "):\n", sep = "") 
		else cat("Values of names which occur in Call:\n") 
		print(v)		
    	}

	options(deparse.cutoff = cutoff)
}
