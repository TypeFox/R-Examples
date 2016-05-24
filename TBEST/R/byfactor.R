byfactor<-function(chrmdata){
	names<-dimnames(chrmdata)[[2]]
	if(length(unique(names))!=ncol(chrmdata))stop("wrong number of variable names")
	dimnames(chrmdata)[[2]]<-sample(names)
	return(chrmdata[,match(names,dimnames(chrmdata)[[2]])])
}
