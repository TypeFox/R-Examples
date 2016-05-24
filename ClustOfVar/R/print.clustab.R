print.clustab <-
function(x, ...){
	if (!inherits(x, "clustab")) 
        	stop("use only with \"clustab\" objects")
	cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  	cat("\n")
	res <- matrix("",2,2)
	colnames(res) <-c("name","description")
	res[1,] <- c("$matCR", "matrix of corrected Rand indices")
	res[2,] <- c("$meanR", "vector of mean corrected Rand indices")
	row.names(res) <- rep("",2) 
	print(res)

}

