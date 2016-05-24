###################################################
### code chunk
###################################################

summary.CHR <- function(object, log.CHR=FALSE, ...) { 
  if(log.CHR==FALSE | log.CHR==F) {
    temp <- list(Call=object$Call, 
                 time=object$time, n.risk=object$n.risk, n.event=object$n.event,
                 CHR1211=object$CHR1211, CHR2111=object$CHR2111, CHR2211=object$CHR2211, 
                 CHR2112=object$CHR2112, CHR2212=object$CHR2212, CHR2221=object$CHR2221,
                 SE1211=object$SE1211, SE2111=object$SE2111, SE2211=object$SE2211,
                 SE2112=object$SE2112, SE2212=object$SE2212, SE2221=object$SE2221) }
  if(log.CHR==TRUE | log.CHR==T) {
    temp <- list(Call=object$Call, 
                 time=object$time, n.risk=object$n.risk, n.event=object$n.event,
                 CHR1211=object$CHR1211.LOG, CHR2111=object$CHR2111.LOG, CHR2211=object$CHR2211.LOG, 
                 CHR2112=object$CHR2112.LOG, CHR2212=object$CHR2212.LOG, CHR2221=object$CHR2221.LOG,
                 SE1211=object$SE1211.LOG, SE2111=object$SE2111.LOG, SE2211=object$SE2211.LOG,
                 SE2112=object$SE2112.LOG, SE2212=object$SE2212.LOG, SE2221=object$SE2221.LOG) }
  class(temp) <- "summary.CHR"
  return(temp)
} 


