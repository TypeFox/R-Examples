summary.clustvar <-
function(object, ...)
{
	x <- object
    	if (!inherits(x, "clustvar")) 
       	stop("use only with \"clustvar\" objects")
	cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
	iter <- x$iter
	if (is.numeric(iter)) cat("number of iterations: ",iter,sep=" ") 
	cat("\n")
	cat("\n")
	n <- x$rec$n
	p <- x$rec$p
  	p1 <- x$rec$p1
  	p2 <- p-p1
	cat("Data:", "\n")
	cat(paste("   number of observations: ",n),sep=" ") 
	cat("\n")
  	if 	((p1!=0)&& (p2==0)) {
  		cat(paste("   number of variables: ",p1),sep=" ")  
  		cat("\n")
	}
  	if 	((p1==0)&& (p2!=0)) {
 		cat(paste("   number of variables: ",p2),sep=" ")  
  		cat("\n")
	}
  	if 	((p1!=0)&& (p2!=0)) {
  		cat(paste("   number of  variables: ",p),sep=" ")
  		cat("\n")
  		cat(paste("        number of numerical variables: ",p1),sep=" ")   
  		cat("\n")
  		cat(paste("        number of categorical variables: ",p2),sep=" ")   
  		cat("\n")
	}
	cat(paste("   number of clusters: ",x$k),sep=" ") 
	cat("\n")
	
	
	for (g in 1:x$k)
	{
		cat("\n")
		cat(paste("Cluster ",g,": "),sep=" ")
		cat("\n")
		print(x$var[[g]],digits=2)
		cat("\n")
	}
	cat("\n")
	cat(paste("Gain in cohesion (in %): ",round(x$E,digits=2)),sep=" ") 
	cat("\n")
}

