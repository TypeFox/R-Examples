`print.subt` <-
function(x,...)
{
cat("\tSubampling t-test",fill=TRUE)
cat("\nOriginal sample sizes:",attr(x,'n1'),'and',
                               attr(x,'n2'),fill=TRUE)
cat("Subsampling method:", if(attr(x,'balanced')) 'balanced' else 'unbalanced',fill=TRUE)
cat("Method to estimate p-value density at 1:",attr(x,'f1method'),fill=TRUE)
cat("Maximum number of subsamples per subsample size configuration:",
                            attr(x,'max.reps'),fill=TRUE)
cat("Total number of subsamples used:", nrow(x),fill=TRUE)
invisible(NULL)
}

