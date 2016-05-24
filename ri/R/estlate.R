estlate <-
function(Y,D,Z,X=NULL,Ypre=NULL,Dpre=NULL,prob=NULL,HT=FALSE) {
	late <- estate(Y,Z,X,Ypre,prob,HT)/estate(D,Z,X,Dpre,prob,HT)
	return(late)
	}
