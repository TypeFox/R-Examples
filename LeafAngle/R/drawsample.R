`drawsample` <-
function(obj, n=25, degrees=FALSE, ...){
	
	mindens <- 0 
	maxdens <- pi/2
	
	if(obj$distribution != "const"){
	RN <- NULL	
	for(i in 1:n){
		OK <- 0
		while(OK < 1) {
			Ts <- runif(1,min = mindens, max = maxdens )
			U <- runif(1,min = 0, max = 1)
						
		if(U*maxdens <= ftheta(Ts,distribution=obj$distribution,distpars=obj$distpars,...) ) {
			OK <- 1
			RN <- c(RN,Ts) 
		}
		}
		
	}
		if(degrees)RN <- RN * 180/pi
	}
	else {
		RN <- rep(obj$distpars[1],n)
	}
	

return(RN)
}

