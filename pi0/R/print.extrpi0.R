`print.extrpi0` <-
function(x,...)
{
cat("\tSubampling-Extrapolation Based Estimation of Proportion of True Null Hypotheses (pi0)",fill=TRUE)
cat("\nOriginal sample sizes:",attr(attr(x,'subt.data'),'n1'),'and',
                               attr(attr(x,'subt.data'),'n2'),fill=TRUE)
cat("Subsampling method:", if(attr(attr(x,'subt.data'),'balanced')) 'balanced' else 'unbalanced',fill=TRUE)
cat("Method to estimate p-value density at 1:",attr(attr(x,'subt.data'),'f1method'),fill=TRUE)
cat("Maximum number of subsamples per subsample size configuration:",
                            attr(attr(x,'subt.data'),'max.reps'),fill=TRUE)
cat("Total number of subsamples used:", nrow(attr(x,'subt.data')),fill=TRUE)
cat("Extrapolation model:", attr(x,'nparm'),'parameters',fill=TRUE)
cat("Nonlinear optimization function:", attr(x,'extrpFUN'),fill=TRUE)
cat("Nonlinear optimization start values:", paste(attr(x,'start.val'),collapse=', '),fill=TRUE)
cat("\nFinal estimate of pi0:", x[1],fill=TRUE)
invisible(NULL)
}

