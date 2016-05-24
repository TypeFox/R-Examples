CrossV <-
function(answered0,subjtheta0,kernel0){
	## This function finds the optimal bandwidth for each item.


	Calc<-function(x,bandwidth){
		blah<-CV(A=bandwidth,B=subjtheta0,C=kerneltog,D=x,E=answered0[,-x])
		## Use Rcpp and C++ to do the heavy lifting.
		blah[which(is.nan(blah) | blah > 1)]<-0
		## Constrain some paramters so algorithm doesn't explode.
		error<-(1-sum(blah*answered0[,x]))**2
		## Return cross validation error.
		return(error)

	}



	CVstat<-function(x){
		bandwidth<-x
		## Do this for a random sample to speed up the algorithm
		CVsamp <- sample(1:length(subjtheta0),round(length(subjtheta0)/10),replace=FALSE)
		return(sum(sapply(CVsamp,Calc,bandwidth=bandwidth)))
	}


		## Pass which kernel we are using to the C++ function.
		if(kernel0=="gaussian"){kerneltog<-1;}
		if(kernel0=="quadratic"){kerneltog<-2;}
		if(kernel0=="uniform"){kerneltog<-3;}

		## Univariate optimization to return the minimum error per item.
		here<-optimize(f=CVstat,interval=c(0,1),tol=.05)
		return(here$minimum)


}

