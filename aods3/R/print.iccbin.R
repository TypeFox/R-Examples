print.iccbin <- function(x, ...){                      
	feat <- x$features
  cat("\nIntra-class correlation")
	cat("\nMethod = ", feat["method"], ", nAGQ = ", feat["nAGQ"], ", M = ", feat["M"],
		", Rho: ", x$rho, "\n\n", sep = "")
  }
