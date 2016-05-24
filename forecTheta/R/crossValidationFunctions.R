
groe <- function(y, forecFunction, g="sAPE", n1=length(y)-10, m=5, H=length(y)-n1, p=1+floor((length(y)-n1)/m), ...){	
	if(n1>=length(y)){ stop("Error in groe function: n1>=length(y)") }
	if(n1<4){ stop("Error in groe function: n1<4") }
	if(m<1){ stop("Error in groe function: m<1") }
	if(H<1){ stop("Error in groe function: H<1") }
	if(p > 1+floor((length(y)-n1)/m)){ stop("ERROR in groe function: p > 1+floor((length(y)-n1)/m)") }
	if(p <= 0){ stop("ERROR in groe function: p <= 0") }
	if(!any( g==c("AE","SE","APE","sAPE") )) stop("Error in lossFunction: this g function has not been implemented.")

	n = length(y)
	
	predictionErrors <- function(i){
		ni = n1+(i-1)*m
		n_pred = min(H,n-ni)
		
		prediction <- forecFunction( as.ts(y[1:ni]), h=n_pred, ... )$mean 	
		
		errors <- errorMetric(obs=y[(ni+1):(ni+n_pred)], forec=prediction, type=g, statistic="N")
    
		if( i < p && n1+i*m < n)
		  errors <- c( errors, predictionErrors(i+1) )
			
		return(errors)
	}

	errors <- predictionErrors(i=1)
	
	return( sum( errors ) )
}


fixOrig <- function(y, forecFunction=stheta, g="sAPE", n1=length(y)-10, ...){
	n = length(y)
	m = n-n1
	H = n-n1
	p=1	
	groe(y=y, forecFunction=forecFunction, g=g, n1=n1, m=m, H=H, p=p, ...)
}


rolOrig <- function(y, forecFunction=stheta, g="sAPE", n1=length(y)-10, ...){
	n = length(y)
	m = 1
	H = n-n1
	p=1+floor((n-n1)/m)	
	groe(y=y, forecFunction=forecFunction, g=g, n1=n1, m=m, H=H, p=p, ...)
}

