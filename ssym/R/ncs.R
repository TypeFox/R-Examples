ncs <-
function(xx, lambda, nknots, all.knots){
	xx <- as.matrix(round(xx,digits=6))
	difv <- as.matrix(as.numeric(levels(factor(xx))))
	if(length(difv) < 3) stop("Variable in natural cubic spline does not have at least three different values!!",call.=FALSE)
	if(!is.numeric(xx)) stop("Variable in natural cubic spline must be numeric!!",call.=FALSE)
	if(missingArg(all.knots))  all.knots <- FALSE

	if(all.knots==FALSE){
		if(missingArg(nknots)){
		  nknots <- floor(length(xx)^(1/3)) + 3
		  xk <- quantile(xx,prob=seq(0,1,length=nknots))
		  difv2 <- as.matrix(as.numeric(levels(factor(xk))))
		  while(length(difv2) < nknots){
		     nknots <- nknots - 1
		     xk <- quantile(xx,prob=seq(0,1,length=nknots))
		     difv2 <- as.matrix(as.numeric(levels(factor(xk))))
		  }
		  if(length(difv) < nknots){
			nknots <- length(difv)
	        xk <- difv
		  }
		}
		else{if(floor(nknots)<3) stop("Number of knots must be an integer >= 3!!",call.=FALSE)
	   		 nknots <- floor(nknots)
   		     xk <- quantile(xx,prob=seq(0,1,length=nknots))
		     difv2 <- as.matrix(as.numeric(levels(factor(xk))))
			 if(length(difv2) < nknots) stop("Too many knots!!",call.=FALSE)
	    }
	}else{nknots <- length(difv)
	      xk <- difv}
	if(!missingArg(lambda)){
	  	if(lambda<=0) stop("Smoothing parameter must be a positive value!!",call.=FALSE)
		else status <- "known"
	}else{status <- "unknown"
	      lambda <- 1
	}

	n <- length(xx)
    m <- nknots
	h <- matrix(0,m-1,1)
	Q <- matrix(0,m,m-2)
	R <- matrix(0,m-2,m-2)
	
	for(i in 1:(m-1)){
	   h[i] <- xk[i+1]-xk[i]
	}
	for(j in 2:(m-1)){
	   Q[j-1,j-1] <- 1/h[j-1]
	   Q[j,j-1] <- -1/h[j-1] -1/h[j]
	   Q[j+1,j-1] <- 1/h[j]
	   R[j-1,j-1] <- (h[j-1]+h[j])/3
	}
	for(j in 2:(m-2)){
	   R[j-1,j] <- (h[j])/6
	   R[j,j-1] <- (h[j])/6
	}
	K <- Q%*%solve(R)%*%t(Q)

	Amin <- matrix(0,n,m)
	Amax <- matrix(0,n,m)	
	Cmin <- matrix(0,n,m)
	Cmax <- matrix(0,n,m)	
	
	for(i in 1:(m-1)){
	   if(i == (m-1)) Ii <- as.vector(xk[i] <= xx & xx <= xk[i+1]) else Ii <- as.vector(xk[i] <= xx & xx < xk[i+1])
	   Amin[,i]   <-  Ii*(xk[i+1] - xx)/h[i]
	   Amax[,i+1] <-  Ii*(xx - xk[i])/h[i]
	   Cmin[,i]   <- -Ii*(xx - xk[i])*(xk[i+1] - xx)*(1 + (xk[i+1] - xx)/h[i])/6
	   Cmax[,i+1] <- -Ii*(xx - xk[i])*(xk[i+1] - xx)*(1 + (xx - xk[i])/h[i])/6	   
	}
	F <- matrix(0,m,m)
	F[2:(m-1),] <- solve(R)%*%t(Q)
	
attr(xx,"K") <- K
attr(xx,"N") <- Amin + Amax + (Cmin + Cmax)%*%F
attr(xx,"status") <- status
attr(xx,"lambda") <- lambda
attr(xx,"knots") <- xk
	Amin <- matrix(0,200,m)
	Amax <- matrix(0,200,m)	
	Cmin <- matrix(0,200,m)
	Cmax <- matrix(0,200,m)	
	xx2 <- seq(min(xx),max(xx),length=200)
	for(i in 1:(m-1)){
	   if(i == (m-1)) Ii <- as.vector(xk[i] <= xx2 & xx2 <= xk[i+1]) else Ii <- as.vector(xk[i] <= xx2 & xx2 < xk[i+1])
	   Amin[,i]   <-  Ii*(xk[i+1] - xx2)/h[i]
	   Amax[,i+1] <-  Ii*(xx2 - xk[i])/h[i]
	   Cmin[,i]   <- -Ii*(xx2 - xk[i])*(xk[i+1] - xx2)*(1 + (xk[i+1] - xx2)/h[i])/6
	   Cmax[,i+1] <- -Ii*(xx2 - xk[i])*(xk[i+1] - xx2)*(1 + (xx2 - xk[i])/h[i])/6	   
	}
	F <- matrix(0,m,m)
	F[2:(m-1),] <- solve(R)%*%t(Q)
attr(xx,"N2") <- Amin + Amax + (Cmin + Cmax)%*%F
xx
}
