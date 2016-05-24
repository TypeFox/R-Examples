`checkX1X2` <-
function(x1,x2,impute=TRUE){
	namex2<-ifelse(impute,"mat.na","newdata")
	if(!is.matrix(x1) | !is.matrix(x2))
		stop("Both data and ",namex2," must be matrix objects.",call.=FALSE)
	if(ncol(x1)!=ncol(x2))
		stop("data and ",namex2," must have the same number of columns.",call.=FALSE)
	m1<-max(x1,na.rm=TRUE)
	m2<-max(x2,na.rm=TRUE)
	if(m1!=m2)
		stop("data and ",namex2," must consist of the same numbers of categories.",call.=FALSE)
	checkCatMat(x1,m1)
	checkCatMat(x2,m2,matname=namex2)
	m1
}

