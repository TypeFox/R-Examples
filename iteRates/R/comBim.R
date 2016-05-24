comBim <-
function(subTax,k){
	#provides a matrix of possible bin sizes given k
	v<-length(subTax)
	comBin<-t(restrictedparts(v,k,include.zero=FALSE))
	return(comBin)
	}

