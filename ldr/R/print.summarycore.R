print.summarycore <-
function(x,...)
{
  	cat("\nCall:\n"); print(x$call);

	ans <- summary.core(x);

	if (identical(x$numdir.test, TRUE))
	{	
		cat("\nEstimated Basis Vectors for Central Subspace:\n");
		
		print(round(x$Gammahat[[x$numdir+1]], digits=4))

		cat("\nInformation Criterion:\n");

		print(ans$IC);

		cat("\nLarge sample likelihood ratio test \n")

		print(ans$LRT); cat("\n");
	}
	else 
	{
		cat("\n\nEstimated Basis Vectors for Central Subspace:\n");

		print(round(x$Gammahat, digits=4)); cat("\n"); 
	}
}
