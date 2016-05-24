print.MultiDEGeneExpr <-
function (x, nlimit=Inf, ...)
{
	cat("Object of class MultiDEGeneExpr (package stringgaussnet)\n\n")
	cat(length(x)," objects of class DEGeneExpr (",paste(names(x),collapse=", "),")\n\n",sep="")
	if (nlimit>length(x)){nlimit<-length(x)}
	for (i in 1:nlimit)
	{
		name<-names(x)[i]
		cat(paste(name,":\n",sep=""))
		print(x[[i]],2)
		cat("\n")
	}
}
