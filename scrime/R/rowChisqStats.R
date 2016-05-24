`rowChisqStats` <-
function(data,cl,compPval=TRUE,asMatrix=TRUE){
	if(!is.matrix(data))
		stop("data has to be a matrix.")
	n.cat<-max(data,na.rm=TRUE)
	checkCatMat(data,n.cat)
	if(missing(cl))
		return(chisqInd(data,n.cat,compPval=compPval,asMatrix=asMatrix))
	else
		return(chisqClass2(data,cl,n.cat,compPval=compPval))
}

