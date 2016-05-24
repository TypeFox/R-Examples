print.variofit <- function (x, digits=max(3L, getOption("digits") - 3L), ...) {
  switch(x$model,
         spherical={name='Spherical model'},
         gaussian={name='Gaussian model'},
         nugget={name='Nugget model'},
         linear={name='Linear model'},
         exponential={name='Exponential model'},
         sill={name='Sill model'},
         periodic={name='Periodic model'},
         hole={name='Hole model'})
           
  if (x$convergence==TRUE)
    conv="algorithm converged"
  else
    conv="algorithm failed to converge"
  
  cat(name)
  cat("\nRoot Mean Square Error (RMSE): ")
  cat(format(x$rmse, digits=digits))
  cat("\nAIC: ")
  cat(format(x$AIC, digits=digits))
  cat(paste("\n\nMaximum likelihood estimates (", conv, "):\n", sep=""))
  print(x$params, digits=digits)
  cat("\nMean bin distance:\n")
  print(x$bins, digits=digits)
  cat("\nObserved values:\n")
  print(x$vario, digits=digits)
  cat("\nPredicted values:\n")
  print(x$fit, digits=digits)  
}
