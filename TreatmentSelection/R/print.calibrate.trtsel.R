print.calibrate.trtsel <-
function( x, ... ) {
cal <- x
HL <- cal$HL.TestStat
p.value <- cal$p.value
Df <- cal$Df

cat("\n")
  cat("  Hosmer - Lemeshow test for model calibration\n")
  cat(" ----------------------------------------------\n\n")

  cat(paste("   Number of Groups:", Df+2, "\n\n"))
  cat("   No Treatment (trt = 0):\n")

  cat("    Test Statistic = "); cat(round(HL[1], 3)); cat(",   DF = "); cat(Df); cat(",   p value = "); cat(p.value[1]); cat("\n\n")

  cat("   Treated (trt = 1):\n")

  cat("    Test Statistic = "); cat(round(HL[2], 3)); cat(",   DF = "); cat(Df); cat(",   p value = "); cat(p.value[2]); cat("\n\n\n")

}
