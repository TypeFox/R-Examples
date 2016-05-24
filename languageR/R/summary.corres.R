`summary.corres` <-
function(object, n = 2, returnList=FALSE, head=TRUE, ...) {
	if (!is(object, "corres")) stop("argument should be a correspondence object")
	if (n > ncol(object@data$origOut$ccor)) stop("n exceeds number of factors")
	#  hier moet tabel met c1 en c2: scores, etc etc., zie p. 113
	#  correlations of columns with factors:  object@data$origOut$ccor[,1:3]
	#  factor projections: object@data$origOut$cproj[,1:3]
  #  contribution: object@data$origOut$ccntr[1:3,]

	cat(paste("\nCall:\ncorres(", object@data$inputName, ")\n\n", sep=""))
    #cat("\nEigenvalues (trivial first eigenvalue removed):\n\n")
	  #if (head) 
    #	cat("    ", head(object@data$eigenvals), " ... \n\n")
	  #else
    #	cat("    ", object@data$eigenvals, "\n\n")
  cat("Eigenvalue rates:\n\n")
	if (head)
        cat("    ", head(object@data$eigenrates/1000), " ... \n\n")
	else
        cat("    ", object@data$eigenrates/1000, "\n\n")

	res =  list(data.frame(coordinates = round(object@data$origOut$cproj[,1],3),
	                       correlations = round(object@data$origOut$ccor[,1],3),
	                       contributions = round(object@data$origOut$ccntr[,1],3)))
	for (i in 2:n) {
		res = c(res, list(  
	            data.frame(coordinates = round(object@data$origOut$cproj[,i],3),
	                       correlations = round(object@data$origOut$ccor[,i],3),
	                       contributions = round(object@data$origOut$ccntr[,i],3)))
		)
	}
	for (i in 1:n) {
		names(res)[i] = paste("Factor ", i, sep = "")
		cat(names(res)[i], "\n\n")
		if (head) {
			print(head(res[[i]]))
			cat("... \n\n")
		}
		else print(res[i])
	}
	invisible(object)
	if (returnList) return(res)
}

