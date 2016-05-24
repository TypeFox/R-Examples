###################################################
### code chunk
###################################################

print.summary.DTR <- function(x, ...) { 
  if (!is.null(x$Call)) { 
    cat("Call: ")
    dput(x$Call)
    cat("\n")
  }
  
  #Combine the results
  temp <- data.frame(time=x$time, n.risk=x$n.risk, n.event=x$n.event,
                     SURV11=x$SURV11, SURV12=x$SURV12, 
                     SURV21=x$SURV21, SURV22=x$SURV22,
                     SE11=x$SE11, SE12=x$SE12, 
                     SE21=x$SE21, SE22=x$SE22)
  print(temp, row.names = FALSE) 
}