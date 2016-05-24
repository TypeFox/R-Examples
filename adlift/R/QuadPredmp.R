"QuadPredmp" <-
function(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g){

#does local quadratic (parabolic) prediction for the point remove based on N 
#points (with intercept as default);
#multiple point consideration

if (length(pointsin)==2){
	Xneighbours<-cbind(X[nbrs],X[nbrs]^2)
}
else{
	Xneighbours<-cbind(X[newnbrs],X[newnbrs]^2)
}

Xremove<-cbind(X[remove],X[remove]^2)

if (intercept){
	Xneighbours<-cbind(1,Xneighbours)
	Xremove<-as.row(c(1,Xremove))
}

if (length(nbrs)>=3){       #possible to have parabola with intercept
	temp<-crossprod(Xneighbours)
	mm<-Rmatsolve(temp)%*%t(Xneighbours)
	coeff1<-NULL

	for (i in 1:length(nbrs)){
		coeff1<-cbind(coeff1,as.row(coefflist[[nbrs[i]]]))
	}
	bhat<-mm%*%matrix(coeff1,ncol=1)
	pred<-Xremove%*%bhat
	weights<-matrix(Xremove,nrow=1)%*%mm	#works out neighbour weights
}

if(length(nbrs)==2){		#if we don't have enough (3) neighbours to define a parabola 
				#automatically, we do a linear prediction
	Xneighbours<-Xneighbours[,1:(1+intercept)]
	Xremove<-as.row(Xremove[,1:(1+intercept)])	
	
	temp<-crossprod(Xneighbours)
	mm<-Rmatsolve(temp)%*%t(Xneighbours)
	coeff1<-NULL
	for (i in 1:length(nbrs)){
		coeff1<-cbind(coeff1,as.row(coefflist[[nbrs[i]]]))
	}
	bhat<-mm%*%matrix(coeff1,ncol=1)
	pred<-Xremove%*%bhat
	weights<-matrix(Xremove,nrow=1)%*%mm	#works out nbr weights
} 

if(length(nbrs)==1){
    mm<-0
    bhat<-1
    weights<-1
    pred<-coefflist[[nbrs]]
}

coeff<-as.row(coeff)

return(list(weights=weights,pred=pred,coeff=coeff))

}
