"LinearPred" <-
function(pointsin,X,coeff,nbrs,remove,intercept,neighbours){

# does local linear prediction for the point remove based on N points (with 
# intercept as default)

Xneighbours<-X[nbrs]
Xneighbours<-as.column(Xneighbours)
Xremove<-X[remove]

if (intercept){
	Xneighbours<-cbind(1,Xneighbours)
	Xremove<-as.row(c(1,Xremove))
}

if (length(nbrs)>=2){
	temp<-crossprod(Xneighbours)
	mm<-Rmatsolve(temp)%*%t(Xneighbours)
	bhat<-mm%*%matrix(coeff[nbrs],ncol=1)
	pred<-Xremove%*%bhat
	weights<-matrix(Xremove,nrow=1)%*%mm  	#works out neighbour weights
}
else{
	mm<-0
	bhat<-1
	weights<-1
	pred<-coeff[nbrs]
}

return(list(weights=weights,pred=pred,coeff=coeff))

}
