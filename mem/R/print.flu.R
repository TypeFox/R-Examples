print.flu <-
function(x, ...){
  threshold<-as.data.frame(round(t(x$pre.post.intervals[1:2,3]),2))
  colnames(threshold)<-c("Pre","Post")
  rownames(threshold)<-"Threshold"
  values<-as.data.frame(round(x$epi.intervals[,4],2))
  colnames(values)<-"Threshold"
  rownames(values)<-paste(c("Medium","High","Very high")," (",round(x$epi.intervals[,1]*100,1),"%)",sep="")
	cat("Call:\n")
	print(x$call)
 	cat("\nEpidemic threshold:\n")
	print(threshold)
 	cat("\nIntensity thresholds:\n")
	print(values)
}
