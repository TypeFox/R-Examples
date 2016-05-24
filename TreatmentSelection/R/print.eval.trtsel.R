print.eval.trtsel <-
function(x, ...){
  
#if we fail to reject the null hypothesis, display a warning
  alpha <- x$test.Null$alpha
  
  if(is.na(x$test.Null$reject)){
    #fitted risks provided...no testing done
    
  }else{
    
  
  test.Null.warning(x$test.Null$reject, x$test.Null$a1a3.pval)

  
  cat("\n\n")
  cat("  Hypothesis test:\n")
  cat(" ------------------\n")
  cat("  H0: No marker-by-treatment interaction")
  cat("\n") 
  cat("                                      ")
  cat(" P value = "); cat(round(x$test.Null$p.value, 5));
  cat("\n")
  cat("                                      ")
  cat(" Z statistic = "); cat(round(x$test.Null$z.statistic, 3)); cat("\n")

  if(!is.null(x$test.Null$a1a3.pval)){
  
  cat("  H0: Treatment effect positivity threshold is outside marker bounds")
  cat("\n") 
  cat("                                      ")
  cat(" P value = "); cat(x$test.Null$a1a3.pval);
  cat("\n")

  }
  }
cat("\n")
 
  cat(paste("  Summary Measure Estimates (with ", round(100*(1-alpha)), "% confidence intervals) \n", sep=""))
  cat(" -----------------------------------------------------------\n")

    if(is.null(x$conf.intervals)) {
  cat("   No confidence intervals calculated...returning NA for lower and upper bounds"); cat("\n\n")
   x$conf.intervals <- matrix( nrow = 2, ncol = 16)
  }

  cat("  Decrease in event rate under marker-based treatment (Theta)\n")
  cat("    Empirical:   ")
  cat(paste(" ", round(x$estimates$Theta.emp,  3), " (",
                 round(unname(x$conf.intervals[1,7]), 3), ",",
                 round(unname(x$conf.intervals[2,7]), 3), ") ", sep = ""))
  cat("\n")
  cat("    Model Based: ")
  cat(paste(" ", round(x$estimates$Theta.mod,  3), " (",
                 round(unname(x$conf.intervals[1,8]), 3), ",",
                 round(unname(x$conf.intervals[2,8]), 3), ") ", sep = ""))
  cat("\n\n")

    cat("  Proportion marker negative:\n")
  cat(paste("   ", round(x$estimates$p.neg,      3), " (",
                   round(unname(x$conf.intervals[1,1]), 3), ",",
                   round(unname(x$conf.intervals[2,1]), 3), ") ", sep = ""))
  cat("\n")
  cat("  Proportion marker positive:\n")
  cat(paste("   ", round(x$estimates$p.pos,      3), " (",
            round(unname(x$conf.intervals[1,2]), 3), ",",
            round(unname(x$conf.intervals[2,2]), 3), ") ", sep = ""))
  cat("\n\n")
  
  cat("  Average benefit of no treatment among marker-negatives (B.neg)\n")
  cat("    Empirical:   ")
  cat(paste(" ", round(x$estimates$B.neg.emp,  3), " (",
                 round(unname(x$conf.intervals[1,3]), 3), ",",
                 round(unname(x$conf.intervals[2,3]), 3), ") ", sep = ""))
  cat("\n")
  cat("    Model Based: ")
  cat(paste(" ", round(x$estimates$B.neg.mod,  3), " (",
                 round(unname(x$conf.intervals[1,4]), 3), ",",
                 round(unname(x$conf.intervals[2,4]), 3), ") ", sep = ""))
  cat("\n\n")


  cat("  Average benefit of treatment among marker-positives (B.pos)\n")
  cat("    Empirical:   ")
  cat(paste(" ", round(x$estimates$B.pos.emp,  3), " (",
                 round(unname(x$conf.intervals[1,5]), 3), ",",
                 round(unname(x$conf.intervals[2,5]), 3), ") ", sep = ""))
  cat("\n")
  cat("    Model Based: ")
  cat(paste(" ", round(x$estimates$B.pos.mod,  3), " (",
                 round(unname(x$conf.intervals[1,6]), 3), ",",
                 round(unname(x$conf.intervals[2,6]), 3), ") ", sep = ""))
  cat("\n\n\n")

  
  
  if(is.null(x$discrete.marker)){
  cat("  Variance in estimated treatment effect: \n  ")
  cat(paste("  ", round(x$estimates$Var.Delta,  3), " (",
                 round(unname(x$conf.intervals[1,9]), 3), ",",
                 round(unname(x$conf.intervals[2,9]), 3), ") ", sep = ""))
  

  cat("\n")
  cat("  Total Gain: \n")
  cat(paste("    ", round(x$estimates$TG,  3), " (",
                 round(unname(x$conf.intervals[1,10]), 3), ",",
                 round(unname(x$conf.intervals[2,10]), 3), ") ", sep = ""))

  cat("\n\n")
  cat("  Marker positivity threshold:  "); cat(round( x$estimates$Marker.Thresh, 3))
  cat("\n\n")
  
  }
  
  cat("  Event Rates:\n")
  cat(" --------------------\n")  
  cat("             Treat all       Treat None    Marker-based Treatment") 
  cat("\n")
  cat(" Empirical: ")
  cat(paste(" ", sprintf("   %.3f      ", round(x$estimates$ER.trt1.emp,  3)),  
            sprintf("     %.3f     ",round(x$estimates$ER.trt0.emp,  3)), 
            sprintf("     %.3f    ", round(x$estimates$ER.mkrbased.emp,  3)), sep=""))
  cat("\n          ")
  cat(paste(" ",sprintf(" (%.3f,%.3f)  ", round(unname(x$c[1,13]), 3), round(unname(x$c[2,13]), 3)), 
            sprintf(" (%.3f,%.3f)  ",     round(unname(x$c[1,11]), 3), round(unname(x$c[2,11]), 3)), 
            sprintf(" (%.3f,%.3f) ",      round(unname(x$c[1,15]), 3),round(unname(x$c[2,15]), 3)), sep = ""))
  
  cat("\n Model Based:   ")
  cat(paste("", sprintf("%.3f      ", round(x$estimates$ER.trt1.mod,  3)),  
            sprintf("     %.3f     ",round(x$estimates$ER.trt0.mod,  3)), 
            sprintf("     %.3f    ", round(x$estimates$ER.mkrbased.mod,  3)), sep=""))
  cat("\n          ")
  cat(paste(" ",sprintf(" (%.3f,%.3f)  ", round(unname(x$c[1,14]), 3), round(unname(x$c[2,14]), 3)), 
            sprintf(" (%.3f,%.3f)  ",     round(unname(x$c[1,12]), 3), round(unname(x$c[2,12]), 3)), 
            sprintf(" (%.3f,%.3f) ",      round(unname(x$c[1,16]), 3),round(unname(x$c[2,16]), 3)), sep = ""))
  
  
  
}
