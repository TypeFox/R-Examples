print.summarypfc <-
function(x,...)
{
  "%^%"<-function(M, pow) 
  { 
    if (prod(dim(M)==list(1,1))) return( as.matrix(M^pow) )
    eigenM = eigen(M) 
    return(eigenM$vectors%*%diag(c(eigenM$values)^pow)%*%t(eigenM$vectors))  
  }
  
  cat("\nCall:\n"); print(x$call)

	ans <- summary.pfc(x)

	if (identical(x$numdir.test, TRUE))
	{	
		cat("\nEstimated Basis Vectors for Central Subspace:\n")

		if ((x$structure == "aniso") | (x$structure=="unstr"))
		{
			print(round(orthonorm(((x$Deltahat[[x$numdir]])%^%(-1))%*%x$Gammahat[[x$numdir]]), digits=4)) 
		} else {
			
			print(round(x$Gammahat[[x$numdir]], digits=4)) 
		}
		
		cat("\nInformation Criterion:\n");
		print(ans$IC);

		cat("\nLarge sample likelihood ratio test \n")
		print(ans$LRT); cat("\n")

		if (is.numeric(x$y)) print(round(ans$Rsq, digits=4)); 	cat("\n")  	
	}
	else 
	{
		cat("\n\nEstimated Basis Vectors for Central Subspace:\n")

		if ((x$structure == "aniso") | (x$structure=="unstr"))
		{
			print(round(orthonorm(((x$Deltahat)%^%(-1))%*%x$Gammahat), digits=4)) 
		} else {
			
			print(round(x$Gammahat, digits=4)); cat("\n") 
		}
	}
}
