`print.cca` <-
function(x, ...){
   cat("\nCanonical Correlation Analysis\n\n")
   #Show the correlations
   cat("Canonical Correlations:\n")
   print(x$corr)
   #Show the coefficients
   cat("\nX Coefficients:\n")
   print(x$xcoef)
   cat("\nY Coefficients:\n")
   print(x$ycoef)
   #Structural correlations
   cat("\nStructural Correlations (Loadings) - X Vars:\n")
   print(x$xstructcorr)
   cat("\nStructural Correlations (Loadings) - Y Vars:\n")
   print(x$ystructcorr)
   cat("\nAggregate Redundancy Coefficients (Total Variance Explained):\n")
   cat("\tX | Y:",x$xrd,"\n")
   cat("\tY | X:",x$yrd,"\n")
   cat("\n")
}

