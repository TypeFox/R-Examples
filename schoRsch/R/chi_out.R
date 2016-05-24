chi_out <- function(chioutput, show.n = FALSE, print = TRUE) {
  # Extract values from chisq.test object
  coefficient <- format(round(chioutput$statistic, 2), nsmall = 2)
  degf <- chioutput$parameter
  p <- format(round(chioutput$p.value, 3), nsmall = 3)
  
  # Whether to show sample size
  if (show.n) {
    n <- paste(", n = ", sum(chioutput$observed),sep="")
  } else {
    n = character(0)
  }
  
  # Format p-value
  p <- paste("p = ", p,sep="")
  p <- gsub("p = 1.000", "p > .999", p, fixed = TRUE)
  p <- gsub("p = 0.000", "p < .001", p, fixed = TRUE)
  p <- gsub("p = 0", "p = ", p, fixed = TRUE)
  
  # Assemble text
  outtext <- paste("chi^2(", degf, n, ") = ", coefficient, ", ", p,sep="")
  
  # Assemble output table
  outtable <- data.frame(Test = paste(chioutput$method, ":",sep=""),
                         Results = outtext)
  
  if (print == TRUE) {
    print(outtable)
  } else {
    outtable
  }
}