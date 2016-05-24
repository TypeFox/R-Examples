###################################################
### code chunk
###################################################

print.summary.CHR <- function(x, ...) { 
  if (!is.null(x$Call)) { 
    cat("Call: ")
    dput(x$Call)
    cat("\n")
  }

  # Combine the results
  temp <- data.frame(
    time=x$time, n.risk=x$n.risk, n.event=x$n.event,
    CHR1211=x$CHR1211, CHR2111=x$CHR2111, CHR2211=x$CHR2211, 
    CHR2112=x$CHR2112, CHR2212=x$CHR2212, CHR2221=x$CHR2221,
    SE1211=x$SE1211, SE2111=x$SE2111, SE2211=x$SE2211,
    SE2112=x$SE2112, SE2212=x$SE2212, SE2221=x$SE2221    
  )
  print(temp, row.names=FALSE)

} 


