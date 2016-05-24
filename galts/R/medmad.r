
medmad.cov <- function(data){
	n <- dim(data)[1]
	m <- dim(data)[2]
	medians <- apply (data, 2, median)
	varcov <- matrix(rep(0,m*m), nrow=m, ncol=m)
	for (i in 1:m){
		for (j in i:m){
			if(i != j){
			    v <- median( (data[,i]-medians[i]) *  (data[,j]-medians[j]) )
				varcov[i,j] <- v
				varcov[j,i] <- v
			}else{
				varcov[i,j] <- median(abs(data[,i]-medians[i]))
			}
		}
	}
	return(varcov)
}

	
medmad <- function(formula, h=NULL, csteps=20){
	x<-model.matrix(formula)
    y<-model.frame(formula)[,1]
    n<-length(y)
    p<-dim(x)[2]
    ind<-rep(0,n)
    if(is.null(h)){
        h<-floor(n/2)+floor((p+1)/2)
    }

   cstep<-function(candidates, csteps){
        cmybetas<-candidates
        indices<-order(abs(y-x%*%cmybetas))[1:p]
        for (i in 1:csteps){
            ols<-lm(y[indices]~x[indices,]-1)
            mybetas<-ols$coefficients
            res<-y-x%*%mybetas
            res2<-abs(res)
            o<-order(res2)
            indices<-sort(o[1:h])
        }
        return(mybetas)
    }


	if( setequal(x[,1] , 1) ) {
		data <- cbind(x[,2:dim(x)[2]],y)
		mdata <- x[,2:dim(x)[2]]
	}else{
		data <- cbind(x,y)
		mdata <- x
	}

	m <- apply (mdata,2,median)
	s<-medmad.cov(mdata)
	dist <- mahalanobis(x=mdata, center=m, cov=s)^2
	ordr <- order(dist)
	w <- rep(0,n)
	w[ordr[1:h]] <- 1
	dta <- as.data.frame(cbind(data, w))
	ols <- lm(formula, data=dta, weights=w)	
	initialbetas <- ols$coefficients
	
	newbetas <- cstep(initialbetas, csteps)
	res <- y-x%*%newbetas
    crit<-sum(sort(res^2)[1:h])

	return( list (coefficients=newbetas, residuals=res, crit=crit))	
	
}



