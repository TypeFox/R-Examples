print.extmod <-
function(x, ...)
			{
		cat("Estimate:\n")
		print(x$Estimate[1])
    	cat("\nlowerCI:\n")
		print(x$lowerCI[1])
		cat("\nupperCI:\n")
		print(x$upperCI[1])
	}
