isb <- function(x,spline.knots=0,knots=NULL,spline.degree=1){
	#DEFINES THE KNOTS (AND REMOVE DUPLICATE VALUES)
	if(is.null(knots)){
		quantiles	<- seq(0, 1, length = spline.knots + 2)
		interval	<- quantile(x, probs = quantiles, names=FALSE)           
		knots		<- unique(interval)
	} 
	
	if(length(knots)<=2){
		#NORMALIZE X IN CASE OF NO INTERIOR KNOTS
		I <- structure(normalize(x,standardize='interval'),
			dim=c(length(x),1))
		return(I)
	} else {
		#COMPUTE SPLINE BASIS
		spline.degree <- min(length(knots)-1,spline.degree)
		I <- structure(isplinebasis(x,knots,spline.degree),
			splineInterval = knots, splineDegree = spline.degree)
		return(I)
	}
}

#=========================================
#=========================================
#=========================================
isplinebasis <- function(x,knots,d){
#Transform a given data into I-splines
	if(is.null(knots) || any(is.na(knots)) || any(diff(knots)==0) || length(knots)<=2)  
	    return(x)
	m	   <- length(knots)
	n	   <- length(x)

	interval <- findInterval(x,knots,all.inside=TRUE)
	M	       <- sapply(sequence(m-1),`==`,interval)
	
	for(i in 2:(d+1)){  
	#CREATE M-SPLINE
		tik	<- c(knots[-1],rep(knots[m],i-2))
		ti	<- c(rep(knots[1],i-2),knots[-m])  
		M	<- M %*% diag(1/(tik-ti))
		
		Dx <- Dt <- array(0,dim=c(m+i-3,m+i-2))
		Dx[1L + 0L:(m+i-4L)*(m+i-2L)] 	<- -1
		Dx[1L:(m+i-3L)*(m+i-2L)] 		    <-  1
		
		Dt[1L + 0L:(m+i-4L)*(m+i-2L)] 	<-  tik
		Dt[1L:(m+i-3L)*(m+i-2L)] 		    <- -ti
				
		M <- (M*x) %*% Dx + M %*% Dt 
	}
	
	#REMOVE INTERCEPT
	M <- M[,-1]
	
	#CREATE I-SPLINE
	S <- array(1,dim=rep(NCOL(M),2))
	S[upper.tri(S)] <- 0
	I <- M%*%S
	
	return(I)
}
