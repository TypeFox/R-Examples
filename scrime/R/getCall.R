`getCall` <-
function(call,n.obs){
	fun<-as.character(call$err.fun)
	dots<-call$"..."
	if(is.null(dots))
		return(paste(fun,"(",n.obs,")",sep=""))
	le<-sapply(dots,length)
	if(any(le>1))
		dots[le>1]<-"..."
	tmp<-paste(names(dots),dots,sep=" = ",collapse=", ")
	paste(fun,"(",n.obs,", ",tmp,")",sep="")
}

