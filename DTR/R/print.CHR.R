###################################################
### code chunk
###################################################

print.CHR <- function(x, ...) {
  if (!is.null(x$Call)) { 
    cat("Call: ")
    dput(x$Call)
    cat("\n")
  }  
  # return the results at a chosen time
  index <- which(abs(x$time-x$time75P)==min(abs(x$time-x$time75P)[which(x$time<=x$time75P)]))
  
  temp <- data.frame(
    comparison=x$comparison,
    time75P=rep(x$time75P, 6),
    CHR=c(x$CHR1211[index], x$CHR2111[index], x$CHR2211[index],
          x$CHR2112[index], x$CHR2212[index], x$CHR2221[index]),
    LCL95=c(x$CHR1211[index]-1.96*x$SE1211[index], x$CHR2111[index]-1.96*x$SE2111[index], 
            x$CHR2211[index]-1.96*x$SE2211[index], x$CHR2112[index]-1.96*x$SE2112[index], 
            x$CHR2212[index]-1.96*x$SE2212[index], x$CHR2221[index]-1.96*x$SE2221[index]),
    UCL95=c(x$CHR1211[index]+1.96*x$SE1211[index], x$CHR2111[index]+1.96*x$SE2111[index], 
            x$CHR2211[index]+1.96*x$SE2211[index], x$CHR2112[index]+1.96*x$SE2112[index], 
            x$CHR2212[index]+1.96*x$SE2212[index], x$CHR2221[index]+1.96*x$SE2221[index]),
    LOGCHR=c(x$CHR1211.LOG[index], x$CHR2111.LOG[index], x$CHR2211.LOG[index],
             x$CHR2112.LOG[index], x$CHR2212.LOG[index], x$CHR2221.LOG[index]),
    LOGLCL95=c(x$CHR1211.LOG[index]-1.96*x$SE1211.LOG[index], x$CHR2111.LOG[index]-1.96*x$SE2111.LOG[index], 
               x$CHR2211.LOG[index]-1.96*x$SE2211.LOG[index], x$CHR2112.LOG[index]-1.96*x$SE2112.LOG[index], 
               x$CHR2212.LOG[index]-1.96*x$SE2212.LOG[index], x$CHR2221.LOG[index]-1.96*x$SE2221.LOG[index]),
    LOGUCL95=c(x$CHR1211.LOG[index]+1.96*x$SE1211.LOG[index], x$CHR2111.LOG[index]+1.96*x$SE2111.LOG[index], 
               x$CHR2211.LOG[index]+1.96*x$SE2211.LOG[index], x$CHR2112.LOG[index]+1.96*x$SE2112.LOG[index], 
               x$CHR2212.LOG[index]+1.96*x$SE2212.LOG[index], x$CHR2221.LOG[index]+1.96*x$SE2221.LOG[index])
  )
  print(temp, row.names = FALSE)
}