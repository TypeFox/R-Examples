print.summarylad <-
function(x,...)
{
  	cat("\nCall:\n"); print(x$call);
	ans <- summary.lad(x);

	if (identical(ans$numdir.test, TRUE))
	{	
		cat("\nEstimated Basis Vectors for Central Subspace:\n");
		print(round(ans$Gammahat[[ans$numdir+1]], digits=4))
		cat("\nInformation Criterion:\n");
		print(ans$IC);
		cat("\nLarge sample likelihood ratio test \n")
		print(ans$LRT); cat("\n");
	}
	else 
	{
		cat("\n\nEstimated Basis Vectors for Central Subspace:\n");
		print(round(ans$Gammahat, digits=4)); cat("\n"); 
	}
}
