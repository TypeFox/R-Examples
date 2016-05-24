summary.Nucleolus <-
function(object, ...){
	#assignInNamespace("summary.Nucleolus", summary.Nucleolus, ns = asNamespace("base"))
   x<-object	
   cat("\n")
   cat("Nucleolus of a", x$mode,"Game","for the given coalitions","\n")
   cat("\n")
   rownames(x$nucleolus)<-NULL
   x$nucleolus<-as.data.frame(x$nucleolus)
   print(x$nucleolus)
}
