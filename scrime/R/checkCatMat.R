`checkCatMat` <-
function(data,n.cat,matname="data"){
	range.data<-range(data,na.rm=TRUE)
	if(range.data[1]!=1 | n.cat!=round(range.data[2]))
		stop(matname," must consist of integers between 1 and ",round(n.cat),".",call.=FALSE)
	if(any(!data[!is.na(data)]%in%1:n.cat))
		stop(matname," must consist of integers between 1 and ",n.cat,".",call.=FALSE)
	if(any(!(1:n.cat)%in%data))
		stop("Some of the values between 1 and ",n.cat," are not in ",matname,".",call.=FALSE)
}

