summary.ShapleyShubik <-
function(object, ...) {
	#assignInNamespace("summary.ShapleyShubik", summary.ShapleyShubik, ns = asNamespace("base"))
   x<-object
   cat("\n")
   cat("Distribution of the agents","\n")
   cat("\n")
   D<-x$Results[1,]
   print(D)
   cat("\n")
   cat("Minimum amount of votes to pass a vote: ",x$Quota,"\n")
   cat("\n")  
   cat("Shapley-Shubik Power Index","\n")
   cat("\n")  
   print(x$Results[3,])
}
