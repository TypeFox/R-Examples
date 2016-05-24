print.hclustvar <-
function(x, ...){
	if (!inherits(x, "hclustvar")) 
        	stop("use only with \"hclustvar\" objects")
	cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  	cat("\n")
  	p <- x$rec$p
  	p1 <- x$rec$p1
  	p2 <- p-p1
  	if 	((p1!=0)&& (p2==0)) {
  		cat(paste("number of variables: ",p1),sep=" ")  
  		cat("\n")
	}
  	if 	((p1==0)&& (p2!=0)) {
 		cat(paste("number of variables: ",p2),sep=" ")  
  		cat("\n")
	}
  	if 	((p1!=0)&& (p2!=0)) {
  		cat(paste("number of  variables: ",p),sep=" ")
  		cat("\n")
  		cat(paste("     number of numerical variables: ",p1),sep=" ")   
  		cat("\n")
  		cat(paste("     number of categorical variables: ",p2),sep=" ")   
  		cat("\n")
	}
	n <- x$rec$n
	cat(paste("number of objects: ",n),sep=" ") 
	cat("\n")
	cat("\n")
	cat("available components: ") 
	res <- c("height","clusmat","merge")
	names(res) <- rep("",3)
	cat(res)
	cat("\n")
}

