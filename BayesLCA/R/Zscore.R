Zscore <-
function(X, fit=NULL, itemprob=NULL, classprob=NULL){
	
	if(is.vector(X)) X<- t(as.matrix(X))
	
	N<- nrow(X)
	M<- ncol(X)
	
	if(is.null(fit)){
		if(is.null(itemprob) | is.null(classprob)){
			stop("Parameters not supplied.")
		} else{
		Tau<- classprob
		Theta<- itemprob			
		}
	} else {
		Tau<- fit$classprob
		Theta<- fit$itemprob
	}
		
	G<- length(Tau)
	
	dum<-array(apply(Theta,1,dbinom, size=1, x=t(X)), dim=c(M,N,G))
	
	Z1<-t(Tau*t(apply(dum, c(2,3), prod)))
	Z<-Z1/rowSums(Z1)
	
	rownames(Z)<- apply(X, 1, paste, collapse="")
	if(any(rowSums(Z)==0))
	{
	  lab1<- rowSums(Z) == 0
	  print(paste("The Z scores for the following data points are undefined with respect to the given parameters:", rownames(Z)[lab1]))
	  }
	Z
	}
