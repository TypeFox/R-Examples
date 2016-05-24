print.fmdsd <-
function(x, mean.print=FALSE, var.print=FALSE, cor.print=FALSE, 
            skewness.print=FALSE, kurtosis.print=FALSE, digits=2, ...)
{
cat("group variable: ",x$group, "\n")
cat("variables: ", x$variables, "\n")
cat("---------------------------------------------------------------\n")
cat("inertia\n"); print(x$inertia, digits=3, ...)
#cat("---------------------------------------------------------------\n")
#cat("contributions\n"); print(x$contributions, ...)
#cat("---------------------------------------------------------------\n")
#cat("qualities\n"); print(x$qualities, ...)
cat("---------------------------------------------------------------\n")
cat("coordinates\n"); print(x$scores, ...)
# Display of means and norms per group (optional)
if (mean.print) 
  {n.group <- length(x$means)
   n.var <- length(x$means[[1]])
   Means <- matrix(unlist(x$means), nrow = n.group, ncol = n.var, byrow = TRUE, 
         dimnames = list(NULL, paste("mean", names(x$means[[1]]), sep = ".")))
   Means <- data.frame(group=names(x$means), Means)
   colnames(Means)[1]=colnames(x$scores)[1]
   st.dev <- unlist(lapply(lapply(x$variances, diag), sqrt))
   Means <- data.frame(Means, matrix(st.dev, nrow = n.group, ncol = n.var, 
         byrow = TRUE, dimnames = list(NULL, paste("sd", names(x$means[[1]])))))
#   Means <- data.frame(Means, norm = x$norm$norm)
   cat("---------------------------------------------------------------\n")
   cat("means, standard deviations and norm by group\n") 
   print(Means, digits=digits, ...)
  }
# Variances per group (optional)
if (var.print) 
  {cat("---------------------------------------------------------------\n")
   cat("variances/covariances by group\n"); print(x$variances, digits=digits, ...)
  }
# Correlation per group (optional)
if (cor.print) 
  {cat("---------------------------------------------------------------\n")
  cat("correlations by group\n"); print(x$correlations, digits=digits, ...)
  }
# Skewness per group (optional)
if (skewness.print) 
  {cat("---------------------------------------------------------------\n")
  cat("skewness coefficients by group\n"); print(x$skewness, digits=digits, ...)
  }
# Kurtosis per group (optional)
if (kurtosis.print) 
  {cat("---------------------------------------------------------------\n")
  cat("kurtosis coefficients by group\n"); print(x$kurtosis, digits=digits, ...)
  }
return(invisible(x))
}
