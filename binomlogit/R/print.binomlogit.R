print.binomlogit <- function(x, ...){
      cat("\nData input:\n\n")
      if(is(x,"binomlogitIndiv")){
        cat("Observations:", x$N)
      } else {
        cat("Observations:", x$t)
      }
      cat("\nCovariates:", x$dims)
      
      cat("\n\n\nMCMC details:\n\n")
      cat("Draws (incl. BI):", x$sim, "draws\n")
      cat("Burn-in:", x$burn, "draws\n\n\n")
      
      cat("Runtime:\n\n")
      cat("Total time: ", x$duration, " sec.\n")
      cat("Time (without BI):", x$duration_wBI, "sec.\n\n")
      
      if(is(x,"binomlogitMH")||is(x,"binomlogitHAM")){
            cat("\nAcceptance rate:", x$rate, "%\n")
      }      
}
