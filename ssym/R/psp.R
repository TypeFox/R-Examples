psp <-
function(xx, lambda, b.order, nknots, diff){
   if(!is.numeric(xx)) stop("Variable in P-spline must be numeric!!",call.=FALSE)
	if(!missingArg(lambda)){
	  	if(lambda<=0) stop("Smoothing parameter must be a positive value!!",call.=FALSE)
		else status <- "known"
	}else{status <- "unknown"
	      lambda <- 1}
   if(missingArg(b.order)) b.order <- 3


   if(missingArg(nknots)){
	   nknots <- min(length(as.matrix(as.numeric(levels(factor(xx))))),floor(length(xx)^(1/3))+3)
	   iknots <- quantile(xx,prob=seq(0,1,length=nknots))
	   difv2 <- as.matrix(as.numeric(levels(factor(iknots))))
	   while(length(difv2) < nknots){
		   nknots <- nknots - 1
		   iknots <- quantile(xx,prob=seq(0,1,length=nknots))
		   difv2 <- as.matrix(as.numeric(levels(factor(iknots))))
	   }
   }else{if(floor(nknots)<=0) stop("Number of internal knots must be a positive integer!!",call.=FALSE)
   		 nknots <- floor(nknots)
		 iknots <- quantile(xx,prob=seq(0,1,length=nknots))
		 difv2 <- as.matrix(as.numeric(levels(factor(iknots))))
		 if(length(difv2) < nknots) stop("Too many knots!!",call.=FALSE)
   }

   iknots <- c(seq(min(xx)*0.9,min(xx)*0.95,length=b.order),iknots,seq(max(xx)*1.05,max(xx)*1.1,length=b.order))

   if(!missingArg(diff)){ 
	  	if(floor(diff)<0) stop("Order of the difference penalty term must be a nonnegative integer!!",call.=FALSE)
		diff <- floor(diff)
   }	
   else diff <- 2

   bspline <- function(x,k,i,m){
	if (m==-1){ res <- as.numeric(x<k[i+1] & x>=k[i])}
	else{z0 <- (x-k[i])/(k[i+m+1]-k[i])
	     z1 <- (k[i+m+2]-x)/(k[i+m+2]-k[i+1])
 	     res <- z0*bspline(x,k,i,m-1)+ z1*bspline(x,k,i+1,m-1)
	}
	res
   }
   total <- nknots + b.order - 1
   B <- matrix(0,length(xx),1)
   for(i in 1:total) B <- cbind(B, bspline(xx,iknots,i,b.order-1))
   xx2 <- B[,-1]
   N <- as.matrix(xx2)
   if(diff > 0) P <- diff(diag(ncol(N)),differences=diff) else P <- diag(ncol(N))
   K <- crossprod(P,P)

   B <- matrix(0,200,1)
   for(i in 1:total) B <- cbind(B, bspline(seq(min(xx),max(xx),length=200),iknots,i,b.order-1))
   xx2 <- B[,-1]
   N2 <- as.matrix(xx2)

x <- as.matrix(xx)   
attr(x,"K") <- K
attr(x,"N") <- N
attr(x,"N2") <- N2
attr(x,"status") <- status
attr(x,"lambda") <- lambda
attr(x,"knots") <- iknots
x
}
