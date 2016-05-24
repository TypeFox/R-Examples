# Produces a summary table for ps object 
print.summary.mniptw <- function(x, ...)
{
      	nSum <- length(x$summaryList)      	
      	for(i in 1:nSum){	
      		cat("Summary for time period ", x$uniqueTimes[i], ": \n")
      		print(x$summaryList[[i]])
      		
      	}

      invisible(x)
     }

