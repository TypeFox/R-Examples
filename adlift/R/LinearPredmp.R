"LinearPredmp" <-
function(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g){

# does local linear prediction for the point remove based on N points (with 
# intercept as default);
# with multiple point consideration

if (length(pointsin)==2){
	Xneighbours<-X[nbrs]
}
else{
	Xneighbours<-X[newnbrs]
}

Xneighbours<-as.column(Xneighbours)
Xremove<-X[remove]

if (intercept){
	Xneighbours<-cbind(1,Xneighbours)
	Xremove<-as.row(c(1,Xremove))
}

if (length(nbrs)>=2){
	temp<-crossprod(Xneighbours)
	mm<-Rmatsolve(temp)%*%t(Xneighbours)
	coeff1<-NULL

	for (i in 1:length(nbrs)){
		coeff1<-cbind(coeff1,as.row(coefflist[[nbrs[i]]]))
	}

	bhat<-mm%*%matrix(coeff1,ncol=1)
	pred<-Xremove%*%bhat
	weights<-matrix(Xremove,nrow=1)%*%mm		#works out neighbour weights
}
else{
	mm<-0
	bhat<-1
	weights<-1
	pred<-coefflist[[nbrs]]
}

coeff<-as.row(coeff)

return(list(weights=weights,pred=pred,coeff=coeff))

}
