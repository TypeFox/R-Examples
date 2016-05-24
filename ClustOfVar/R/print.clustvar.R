print.clustvar <-
function(x, ...){
	if (!inherits(x, "clustvar")) 
        	stop("use only with \"clustvar\" objects")
	cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
	iter <- x$iter
	if (is.numeric(iter)) cat("number of iterations: ",iter,sep=" ") 
  	cat("\n")
 	cat("\n")
	res <- matrix("",7,2)
	colnames(res) <-c("name","description")
	res[1,] <- c("$var", "list of variables in each cluster")
	res[2,] <- c("$sim", "similarity matrix in each cluster")
	res[3,] <- c("$cluster", "cluster memberships")
	res[4,] <- c("$wss", "within-cluster sum of squares")
	res[5,] <- c("$E", "gain in cohesion (in %)")
	res[6,] <- c("$size", "size of each cluster")
	res[7,] <- c("$scores", "score of each cluster")
	row.names(res) <- rep("",7)

print(res)
}

