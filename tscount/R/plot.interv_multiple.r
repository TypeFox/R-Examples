plot.interv_multiple <- function(x, ...){
  timser <- x$fit_H0$ts
  intervention <- x$interventions
  plot(timser, type="n", main="Time series with detected interventions", xlab="Time", ylab="Value", ...)
  if(nrow(intervention)==0){
    lines(timser)
    legend("topleft", legend="No interventions detected    ", bg="white")
  }else{  
    lines(x=time(timser)[seq(along=timser)], y=x$ts_cleaned[[nrow(intervention)]], lty="dashed", col="blue")
    abline(v=time(timser)[intervention[, "tau"]], col="red")
    lines(timser)
    legend("topleft", legend=paste("Intervention ", rownames(intervention), ": tau=", intervention[,"tau"], ", delta=", intervention[,"delta"], ", size=", round(intervention[,"size"],2), "    ", sep=""), bg="white")
  }
}
