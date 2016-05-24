"QuadPred" <-
function(pointsin,X,coeff,nbrs,remove,intercept,neighbours){

#does local quadratic (parabolic) prediction for the point remove based on N 
#points (with intercept as default);

Xneighbours<-cbind(X[nbrs],X[nbrs]^2)
Xremove<-cbind(X[remove],X[remove]^2)

if (intercept){
	Xneighbours<-cbind(1,Xneighbours);
	Xremove<-as.row(c(1,Xremove));
}

if (length(nbrs)>=3){       #possible to have parabola with intercept

	temp<-crossprod(Xneighbours)
	mm<-Rmatsolve(temp)%*%t(Xneighbours)
	bhat<-mm%*%matrix(coeff[nbrs],ncol=1)
	pred<-Xremove%*%bhat
	weights<-matrix(Xremove,nrow=1)%*%mm  #works out neighbour weights
}


if(length(nbrs)==2){
	
	#if we don't have enough (3) neighbours to define a parabola 
	#automatically, we do a linear prediction
	
	Xneighbours<-Xneighbours[,1:(1+intercept)]
	Xremove<-as.row(Xremove[,1:(1+intercept)])	

	temp<-crossprod(Xneighbours)
	mm<-Rmatsolve(temp)%*%t(Xneighbours)

	bhat<-mm%*%matrix(coeff[nbrs],ncol=1)
	
	pred<-Xremove%*%bhat
		
	weights<-matrix(Xremove,nrow=1)%*%mm  #works out nbr weights
				
}

if(length(nbrs)==1){
	mm<-0
	bhat<-1
	weights<-1
	pred<-coeff[nbrs]
}

return(list(weights=weights,pred=pred,coeff=coeff))

}
