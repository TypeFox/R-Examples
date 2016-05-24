print.etasclass<-function(x,...){
	  if(class(x)!="etasclass")stop(" argument x must be an etasclass object")

cat("Call:","\n","\n")
print(x$this.call)
cat("\n","\n","\n")
cat(x$description,"\n")
cat("Number of observations            ",length(x$cat$time),"\n")

cat("ETAS Parameters:","\n")
ris=x$params
print(round(ris,6))
}
