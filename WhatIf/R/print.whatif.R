print.whatif <- function(x, print.dist = FALSE, print.freq = FALSE, ...)  {

   x$in.hull <- as.character(x$in.hull)

    #LOCAL VARIABLES
    m <- length(x$in.hull)
 
    #PRINT ON SCREEN
    cat("\nCall:  ", deparse(x$call), sep="\n")
    cat("\n")
    cat("Counterfactual in Convex Hull, True or False:\n")
    prmatrix(data.frame(Counterfactual = seq(1, m, by = 1), "In Hull" = x$in.hull, 
      check.names = FALSE), rowlab = rep("", m), quote = FALSE, right = TRUE)
    cat("\n")
    cat("Percent Data Nearby Counterfactual:\n")
    prmatrix(data.frame(Counterfactual = seq(1, m, by = 1), "Percent Nearby" = x$sum.stat, 
      check.names = FALSE), rowlab = rep("", m), quote = FALSE, right = TRUE)
    cat("\n")
    cat("Geometric Variance of Covariates:  ", x$geom.var, sep = "\n")
    cat("\n")
    if (print.dist)  {
      if (is.null(x$dist))  {
        print("No distance matrix returned to print")
      }  else  {
        cat("Distances of Counterfactual to Data Points:\n")
        prmatrix(cbind(Counterfactual = seq(1, m, by = 1), x$dist), rowlab = 
          rep("", m))
        cat("\n")
      }
    }
    if (print.freq)  {
      cat("Cumulative Frequencies of Distances:\n")
      prmatrix(cbind(Counterfactual = seq(1, m, by = 1), x$cum.freq), rowlab = 
        rep("", m))
      cat("\n")
   }
   
    return(invisible(x))
}
