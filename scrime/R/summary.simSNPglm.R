`summary.simSNPglm` <-
function(object,digits=3,...){
	print(object)
	if(!is.null(object$err))
		return(invisible())
	tab<-table(object$prob,object$y)
	out<-data.frame(Probability=round(as.numeric(rownames(tab)),digits),
		"    1"=tab[,"1"],"    0"=tab[,"0"],check.names=FALSE)
	rownames(out)<-0:(nrow(out)-1)
	cat("\n","Number of Cases (1) and Controls (0):\n",sep="")
	print(out)
	invisible(out)
}

