summary.ShapleyValue <-
function(object, ...) {
	#assignInNamespace("summary.ShapleyValue", summary.ShapleyValue, ns = asNamespace("base"))
	x<-object
   cat("\n")
   cat("Shapley Value for the given game","\n")
   cat("\n")
   print(as.data.frame(x$SV))
}
