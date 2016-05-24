print.MultiNetworks <-
function (x, ...)
{
	cat("Object of class MultiNetworks (package stringgaussnet)\n\n")
	cat(length(x)," object(s) of class DEGeneExpr used (",paste(names(x),collapse=", "),")\n\n",sep="")
	UseMethods<-grep("(^STRING$)|(^SIMoNe$)|(^WGCNA$)",names(x[[1]]),value=T)
	cat(length(UseMethods)," method(s) of network creation used (",paste(UseMethods,collapse=", "),")\n\n",sep="")
	FactorNets<-grep("FactorNetworks",names(x[[1]]),value=T)
	if (length(FactorNets)==0){cat("No factor entered by the user\n")}
	if (length(FactorNets)>0)
	{
		Levels<-names(x[[1]][[FactorNets[1]]])
		cat("A factor with ",length(Levels)," levels has been entered by the user (",paste(Levels,collapse=", "),")\n",sep="")
	}
}
