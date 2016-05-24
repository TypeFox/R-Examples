print.compare.trtsel <-
function(x, ...){
  
  confperc <- round(100*(1-x$alpha))

  cat("                      Summary Measure Estimates \n")
  cat(paste("                    (with ", confperc ,"% confidence intervals) \n"))

  cat("\n               marker 1    |    marker 2    |   difference    (p-value)\n" )
  cat(" ------------------------------------------------------------------------\n\n")


 

  cat("Decrease in event rate under marker-based treatment (Theta)\n")
   bootN <- x$bootstraps  

   if(bootN==0){ 
     x$p.values = rep(NA, 10)
     cat("   No confidence intervals calculated...returning NA for lower and upper bounds"); cat("\n\n")
     x$ci.marker1 <- matrix( nrow = 2, ncol = 10)
     x$ci.marker2 <- matrix( nrow = 2, ncol = 10)
     x$ci.diff <- matrix( nrow = 2, ncol = 1)
   }

     tmp.pval <- x$p.values[7]
  if(is.element(tmp.pval, 0) & bootN>0) tmp.pval <- paste("<", 1/bootN)



  cat(" Empirical: ")
  cat(paste(" ", sprintf("   %.3f      |", round(x$estimates.marker1$Theta.emp,  3)),  
                 sprintf("     %.3f     |",round(x$estimates.marker2$Theta.emp,  3)), 
                 sprintf("     %.3f    ", round(x$estimates.diff$Theta.emp,  3)), 
                 "     (", tmp.pval,")", sep=""))


                 
  cat("\n          ")
  cat(paste(" ",sprintf(" (%.3f,%.3f) |", round(unname(x$ci.marker1[1,7]), 3), round(unname(x$ci.marker1[2,6]), 3)), 
                sprintf(" (%.3f,%.3f) |", round(unname(x$ci.marker2[1,7]), 3), round(unname(x$ci.marker2[2,6]), 3)), 
                sprintf(" (%.3f,%.3f) ", round(unname(x$ci.diff[1,7]), 3),round(unname(x$ci.diff[2,6]	), 3)), sep = ""))
     tmp.pval <- x$p.values[8]
  if(is.element(tmp.pval, 0)) tmp.pval <- paste("<", 1/bootN)
  
  cat("\n")
  cat(" Model Based:")  
  cat(paste(" ", sprintf("  %.3f      |", round(x$estimates.marker1$Theta.mod,  3)),  
                 sprintf("     %.3f      |",round(x$estimates.marker2$Theta.mod,  3)), 
                 sprintf("     %.3f    ", round(x$estimates.diff$Theta.mod,  3)), 
                 "     (", tmp.pval,")", sep=""))


                 
  cat("\n          ")
  cat(paste(" ",sprintf(" (%.3f,%.3f)  |", round(unname(x$ci.marker1[1,8]), 3), round(unname(x$ci.marker1[2,7]), 3)), 
                sprintf(" (%.3f,%.3f)  |", round(unname(x$ci.marker2[1,8]), 3), round(unname(x$ci.marker2[2,7]), 3)), 
                sprintf(" (%.3f,%.3f) ", round(unname(x$ci.diff[1,8]), 3),round(unname(x$ci.diff[2,7]), 3)), sep = ""))

  cat("\n")
  cat("\n")

  cat("Proportion marker negative:\n")
     tmp.pval <- x$p.values[1]
  if(is.element(tmp.pval, 0)) tmp.pval <- paste("<", 1/bootN)
  
  cat("            ")  
  cat(paste(" ", sprintf("   %.3f      |", round(x$estimates.marker1$p.neg,  3)),  
                 sprintf("     %.3f      |",round(x$estimates.marker2$p.neg,  3)), 
                 sprintf("     %.3f    ", round(x$estimates.diff$p.neg,  3)), 
                 "     (", tmp.pval,")", sep=""))


                 
  cat("\n          ")
  cat(paste(" ",sprintf(" (%.3f,%.3f)  |", round(unname(x$ci.marker1[1,1]), 3), round(unname(x$ci.marker1[2,1]), 3)), 
                sprintf(" (%.3f,%.3f)  |", round(unname(x$ci.marker2[1,1]), 3), round(unname(x$ci.marker2[2,1]), 3)), 
                sprintf(" (%.3f,%.3f) ", round(unname(x$ci.diff[1,1]), 3),round(unname(x$ci.diff[2,1]), 3)), sep = ""))
  cat("\n")

  
  cat("Proportion marker positive:\n")
  tmp.pval <- x$p.values[2]
  if(is.element(tmp.pval, 0)) tmp.pval <- paste("<", 1/bootN)
  
  cat("            ")  
  cat(paste(" ", sprintf("   %.3f      |", round(x$estimates.marker1$p.pos,  3)),  
            sprintf("     %.3f      |",round(x$estimates.marker2$p.pos,  3)), 
            sprintf("     %.3f    ", round(x$estimates.diff$p.pos,  3)), 
            "     (", tmp.pval,")", sep=""))
  
  
  
  cat("\n          ")
  cat(paste(" ",sprintf(" (%.3f,%.3f)  |", round(unname(x$ci.marker1[1,2]), 3), round(unname(x$ci.marker1[2,2]), 3)), 
            sprintf(" (%.3f,%.3f)  |", round(unname(x$ci.marker2[1,2]), 3), round(unname(x$ci.marker2[2,2]), 3)), 
            sprintf(" (%.3f,%.3f) ", round(unname(x$ci.diff[1,2]), 3),round(unname(x$ci.diff[2,2]), 3)), sep = ""))
  
  cat("\n")
  cat("\n")

     tmp.pval <- x$p.values[3]
  if(is.element(tmp.pval, 0)) tmp.pval <- paste("<", 1/bootN)
  
  cat("Average benefit of no treatment among marker-negatives (B.neg)\n")
  cat(" Empirical: ")
  cat(paste(" ", sprintf("   %.3f      |", round(x$estimates.marker1$B.neg.emp,  3)),  
                 sprintf("     %.3f     |",round(x$estimates.marker2$B.neg.emp,  3)), 
                 sprintf("     %.3f    ", round(x$estimates.diff$B.neg.emp,  3)), 
                 "     (", tmp.pval,")", sep=""))


                 
  cat("\n          ")
  cat(paste(" ",sprintf(" (%.3f,%.3f) |", round(unname(x$ci.marker1[1,3]), 3), round(unname(x$ci.marker1[2,3]), 3)), 
                sprintf(" (%.3f,%.3f) |", round(unname(x$ci.marker2[1,3]), 3), round(unname(x$ci.marker2[2,3]), 3)), 
                sprintf(" (%.3f,%.3f) ", round(unname(x$ci.diff[1,3]), 3),round(unname(x$ci.diff[2,3]), 3)), sep = ""))
  tmp.pval <- x$p.values[4]
  if(is.element(tmp.pval, 0)) tmp.pval <- paste("<", 1/bootN)
  
  cat("\n")
  cat(" Model Based:")  
  cat(paste(" ", sprintf("  %.3f      |", round(x$estimates.marker1$B.neg.mod,  3)),  
                 sprintf("     %.3f      |",round(x$estimates.marker2$B.neg.mod,  3)), 
                 sprintf("     %.3f    ", round(x$estimates.diff$B.neg.mod,  3)), 
                 "     (", tmp.pval,")", sep=""))


                 
  cat("\n          ")
  cat(paste(" ",sprintf(" (%.3f,%.3f)  |", round(unname(x$ci.marker1[1,4]), 3), round(unname(x$ci.marker1[2,4]), 3)), 
                sprintf(" (%.3f,%.3f)  |", round(unname(x$ci.marker2[1,4]), 3), round(unname(x$ci.marker2[2,4]), 3)), 
                sprintf(" (%.3f,%.3f) ", round(unname(x$ci.diff[1,4]), 3),round(unname(x$ci.diff[2,4]), 3)), sep = ""))

  cat("\n")
  cat("\n")

    tmp.pval <- x$p.values[5]
  if(is.element(tmp.pval, 0)) tmp.pval <- paste("<", 1/bootN)
  

  cat("Average benefit of treatment among marker-positives (B.pos)\n")
  cat(" Empirical: ")
  cat(paste(" ", sprintf("   %.3f      |", round(x$estimates.marker1$B.pos.emp,  3)),  
                 sprintf("     %.3f     |",round(x$estimates.marker2$B.pos.emp,  3)), 
                 sprintf("     %.3f    ", round(x$estimates.diff$B.pos.emp,  3)), 
                 "     (", tmp.pval,")", sep=""))


                 
  cat("\n          ")
  cat(paste(" ",sprintf(" (%.3f,%.3f) |", round(unname(x$ci.marker1[1,5]), 3), round(unname(x$ci.marker1[2,5]), 3)), 
                sprintf(" (%.3f,%.3f) |", round(unname(x$ci.marker2[1,5]), 3), round(unname(x$ci.marker2[2,5]), 3)), 
                sprintf(" (%.3f,%.3f) ", round(unname(x$ci.diff[1,5]), 3),round(unname(x$ci.diff[2,5]), 3)), sep = ""))
   tmp.pval <- x$p.values[6]
  if(is.element(tmp.pval, 0)) tmp.pval <- paste("<", 1/bootN)
  
  cat("\n")
  cat(" Model Based:")  
  cat(paste(" ", sprintf("  %.3f      |", round(x$estimates.marker1$B.pos.mod,  3)),  
                 sprintf("     %.3f      |",round(x$estimates.marker2$B.pos.mod,  3)), 
                 sprintf("     %.3f    ", round(x$estimates.diff$B.pos.mod,  3)), 
                 "     (",tmp.pval,")", sep=""))


                 
  cat("\n          ")
  cat(paste(" ",sprintf(" (%.3f,%.3f)  |", round(unname(x$ci.marker1[1,6]), 3), round(unname(x$ci.marker1[2,6]), 3)), 
                sprintf(" (%.3f,%.3f)  |", round(unname(x$ci.marker2[1,6]), 3), round(unname(x$ci.marker2[2,6]), 3)), 
                sprintf(" (%.3f,%.3f) ", round(unname(x$ci.diff[1,6]), 3),round(unname(x$ci.diff[2,6]), 3)), sep = ""))

  cat("\n\n\n")
  tmp.pval <- x$p.values[9]
  if(is.element(tmp.pval, 0)) tmp.pval <- paste("<", 1/bootN)
  
    cat("Variance in estimated treatment effect : \n  ")
  cat("          ")  
  cat(paste(" ", sprintf("   %.3f      |", round(x$estimates.marker1$Var.Delta,  3)),  
                 sprintf("     %.3f      |",round(x$estimates.marker2$Var.Delta,  3)), 
                 sprintf("     %.3f    ", round(x$estimates.diff$Var.Delta,  3)), 
                 "     (", tmp.pval,")", sep=""))  
  cat("\n          ")
  cat(paste(" ",sprintf(" (%.3f,%.3f)  |", round(unname(x$ci.marker1[1,9]), 3), round(unname(x$ci.marker1[2,9]), 3)), 
                sprintf(" (%.3f,%.3f)  |", round(unname(x$ci.marker2[1,9]), 3), round(unname(x$ci.marker2[2,9]), 3)), 
                sprintf(" (%.3f,%.3f) ", round(unname(x$ci.diff[1,9]), 3),round(unname(x$ci.diff[2,9]), 3)), sep = ""))

  cat("\n\n")

  tmp.pval <- x$p.values[10]
  if(is.element(tmp.pval, 0)) tmp.pval <- paste("<", 1/bootN)
  
  cat("Total Gain: \n")
  cat("          ")  
  cat(paste(" ", sprintf("     %.3f      |", round(x$estimates.marker1$TG,  3)),  
                 sprintf("     %.3f      |",round(x$estimates.marker2$TG,  3)), 
                 sprintf("     %.3f    ", round(x$estimates.diff$TG,  3)), 
                 "     (", tmp.pval,")", sep=""))  
  cat("\n          ")
  cat(paste(" ",sprintf(" (%.3f,%.3f)  |", round(unname(x$ci.marker1[1,10]), 3), round(unname(x$ci.marker1[2,10]), 3)), 
                sprintf(" (%.3f,%.3f)  |", round(unname(x$ci.marker2[1,10]), 3), round(unname(x$ci.marker2[2,10]), 3)), 
                sprintf(" (%.3f,%.3f) ", round(unname(x$ci.diff[1,10]), 3),round(unname(x$ci.diff[2,10]), 3)), sep = ""))

  cat("\n\n")  
}
