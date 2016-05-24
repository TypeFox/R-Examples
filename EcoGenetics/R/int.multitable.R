#' Table construction for multiple test result.
#' 
#' @param result Results of \code{\link{int.random.test}}.
#' @param test Test class of the parameter "result" (permutation or bootstrap).
##' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @keywords internal

int.multitable <- function(result, 
													 test = c("permutation", "bootstrap")) {
	
	test <- match.arg(test)
	
	if(test == "permutation") {
		
		res <- data.frame(matrix(nrow = length(result), ncol = 4))
		res [, 1] <- sapply(result, "[[", 3)
		res [, 2] <- sapply(result, "[[", 4)
		res [, 3] <- sapply(result, "[[", 2)
		res [, 4] <- sapply(result, "[[", 6)

		colnames(res) <- c("obs", "exp", "alter", "pval")
		rownames(res) <- names(result)
		res<-list(analysis = result[[1]]$analysis, 
							nsim = result[[1]]$nsim,  
							results = res)
		res
		
	} else if(test == "bootstrap") {
		
		res <- data.frame(matrix(nrow = length(result), ncol = 3))
		res [, 1] <- sapply(result, "[[", 2)
		temp <- lapply(result, "[[", 4)
		temp <- sapply(temp, c)
		res [, 2] <- temp[1, ]
		res [, 3] <- temp[2, ]
		colnames(res) <- c("obs", "lwr", "uppr")
		rownames(res) <- names(result)
		res<-list(analysis = result[[1]]$analysis, 
							nsim = result[[1]]$nsim,  
							results = res)
		res
		
	}
	
	res
	
}
