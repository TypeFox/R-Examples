summary.etasclass<-function(object,...){
x=object
cat("Call:","\n","\n")
print(x$this.call)
cat("\n","\n","\n")
cat(x$description,"\n")
cat("Execution started:                 ",format(x$time.start),"\n")
cat("Elapsed time of execution (hours)  ",as.numeric(x$time.elapsed,units="hours"),"\n")
cat("Number of observations            ",length(x$cat$time),"\n")
cat("Magnitude threshold               ",x$magn.threshold,"\n")

## output also declustering=TRUE|FALSE. the following is only if declustering=TRUE
#declustering=TRUE
cat("declustering                       ",x$declustering,"\n")
if (x$declustering){
cat("Number of declustering iterations  ",x$iter,"\n")
cat("Kind of declustering               ",ifelse(x$thinning,"thinning","weighting"),"\n")
cat("flp                                ",x$flp,"\n")

cat("sequence of AIC values for each iteration","\n")
cat(x$AIC.iter,"\n","\n")
}
cat("-------------------------------------------------------","\n","\n")

cat("ETAS Parameters:","\n")
ris=cbind(x$params,x$sqm)
colnames(ris)=c("      Estimates","      std.err.")
print(round(ris,6))
cat("-------------------------------------------------------","\n")
}
