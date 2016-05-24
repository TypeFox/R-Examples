MultiDEGeneExpr.default <-
function (...)
{
	optionValue<-getOption("stringsAsFactors")
	options(stringsAsFactors=F)
	MultiDEGeneExpr <- list(...)
	for (i in 1:length(MultiDEGeneExpr))
	{
		if (class(MultiDEGeneExpr[[i]])!="DEGeneExpr")
		{
			stop("Bad classes used")
		}
	}
	if(is.null(names(MultiDEGeneExpr))){names(MultiDEGeneExpr)<-paste("List",1:length(MultiDEGeneExpr),sep="")}
	class(MultiDEGeneExpr)<-"MultiDEGeneExpr"
	options(stringsAsFactors=optionValue)
	return(MultiDEGeneExpr)
}
