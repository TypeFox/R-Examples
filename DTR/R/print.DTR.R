###################################################
### code chunk
###################################################

print.DTR <- function(x, ...) {
  if (!is.null(x$Call)) { 
    cat("Call: ")
    dput(x$Call)
    cat("\n")
  }  
  # return the median survival time and the 95% Confidence interval
  temp <- data.frame(
    DTR=x$DTR,
    records=x$records,
    events=x$events,
    median=c(min(x$time[which(abs(x$SURV11-0.5)==min(abs(x$SURV11-0.5)[which(x$SURV11<=0.5)]))]),
             min(x$time[which(abs(x$SURV12-0.5)==min(abs(x$SURV12-0.5)[which(x$SURV12<=0.5)]))]),
             min(x$time[which(abs(x$SURV21-0.5)==min(abs(x$SURV21-0.5)[which(x$SURV21<=0.5)]))]),
             min(x$time[which(abs(x$SURV22-0.5)==min(abs(x$SURV22-0.5)[which(x$SURV22<=0.5)]))])),
    LCL95=c(min(x$time[which(abs(x$SURV11-1.96*x$SE11-0.5)==min(abs(x$SURV11-1.96*x$SE11-0.5)[which((x$SURV11-1.96*x$SE11)<=0.5)]))]),
            min(x$time[which(abs(x$SURV12-1.96*x$SE12-0.5)==min(abs(x$SURV12-1.96*x$SE12-0.5)[which((x$SURV12-1.96*x$SE12)<=0.5)]))]),
            min(x$time[which(abs(x$SURV21-1.96*x$SE21-0.5)==min(abs(x$SURV21-1.96*x$SE21-0.5)[which((x$SURV21-1.96*x$SE21)<=0.5)]))]),
            min(x$time[which(abs(x$SURV22-1.96*x$SE22-0.5)==min(abs(x$SURV22-1.96*x$SE22-0.5)[which((x$SURV22-1.96*x$SE22)<=0.5)]))])),
    UCL95=c(min(x$time[which(abs(x$SURV11+1.96*x$SE11-0.5)==min(abs(x$SURV11+1.96*x$SE11-0.5)[which((x$SURV11+1.96*x$SE11)<=0.5)]))]),
            min(x$time[which(abs(x$SURV12+1.96*x$SE12-0.5)==min(abs(x$SURV12+1.96*x$SE12-0.5)[which((x$SURV12+1.96*x$SE12)<=0.5)]))]),
            min(x$time[which(abs(x$SURV21+1.96*x$SE21-0.5)==min(abs(x$SURV21+1.96*x$SE21-0.5)[which((x$SURV21+1.96*x$SE21)<=0.5)]))]),
            min(x$time[which(abs(x$SURV22+1.96*x$SE22-0.5)==min(abs(x$SURV22+1.96*x$SE22-0.5)[which((x$SURV22+1.96*x$SE22)<=0.5)]))]))
  )
  print(temp, row.names = FALSE)
}