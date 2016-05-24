## calculates spatial correlation
# input:
# X: dataset as nxp matrix
# location: either a numeric with the location or a character of the used location estimator
# 	    possible are: 2dim-median, 1dim-median, pdim-median, mean
# scale: either a numeric with the scale or a character of the used scale estimator
#	 possible are: mad, Qn, sd
# standardized: should the data be standardized first, so that marginal variances are equal?
# output: estimated spatial correlation matrix 

sscor <- function(X,location=c("2dim-median","1dim-median","pdim-median","mean"),scale=c("mad","Qn","sd"),standardized=TRUE,pdim=FALSE,...) {

# protective measures

if(is.character(location)){
	if(length(location)>1){
		if(pdim){location <- "pdim-median"} else {
			location <- location[1]}
		}
	if(sum(location==c("2dim-median","1dim-median","pdim-median","mean"))==0) {
		warning("This location estimator is not implemented. Two dimensional spatial median is used instead.")
		location <- "2dim-median"
		}
	if(pdim & location=="2dim-median") {
		warning("The two dimensional spatial median is not appropriate for the p dimensional spatial sign correleation. P dimensional median is used instead.")
		location <- "pdim-median"
		}
	}
if(is.character(scale)){
	if(length(scale)>1){
		scale <- scale[1]
		}
	if(sum(scale==c("mad","Qn","sd"))==0) {
		warning("This scale estimator is not implemented. MAD is used instead.")
		scale <- "mad"
		}
	}

p <- length(X[1,])
n <- length(X[,1])



if (!pdim) {
	schlussel <- combn(1:p,2)
	nn <- length(schlussel[1,])	# indices for all bivariate correlations
	}

# calculation of marginal variances

if (standardized==FALSE)
	{scale.est <- rep(1,p)}	else{
	if(is.numeric(scale)){
		if (length(scale)==p)
			{scale.est <- scale}	else{
			warning("Scale vector has false dimension. Scale will be estimated using MAD")
			scale <- "mad"
			}
		}	
	if(is.character(scale)) {
			if(scale=="mad") 
				{scale.est <- apply(X,2,mad)}
			if(scale=="Qn")
				{scale.est <- apply(X,2,Qn)}
			if(scale=="sd")
				{scale.est <- apply(X,2,sd)}
		}			
	}				
if (sum(scale.est==0)>0) {
	warning("Some marginal variances are 0. Standardization is not possible.")
	scale.est <- rep(1,p)
	}

# standardisation

X <- t(t(X)/scale.est)

# location estimation

if(is.numeric(location)) {
	if(length(location)==p)
		{location.est <- location} else
			{
			warning("Location vector has false dimension. Location will be estimated using two dimensional spatial median.")
			location <- "2dim-median"
			}
	}
if(is.character(location)){
		if (location=="pdim-median")
			{
			location.est <- try(l1median_BFGS(X)$par,silent=TRUE)
				if (inherits(location.est,"try-error")){
					warning("Calculation of spatial median failed, marginal median will be used instead.")
					location <- "1dim-median"
				} 
			}
		if (location=="mean")
			{location.est <- apply(X,2,mean)}
		if (location=="2dim-median")
			{	
				location.est <- matrix(ncol=2,nrow=nn)
				for(i in 1:nn) {
					location.est[i,] <- try(l1median_BFGS(cbind(X[,schlussel[1,i]],X[,schlussel[2,i]]))$par,silent=TRUE)
					if (inherits(location.est,"try-error")){
						warning("Calculation of spatial median failed. Therefore sometimes marginal median will be used.")
						location.est[i,] <- c(median(X[,schlussel[1,i]]),median(X[,schlussel[2,i]]))
					} 
				}
			}
		if (location=="1dim-median")
			{location.est <- apply(X,2,median)}
		}

if (!pdim) {
	if (length(location.est)==p) {
		location2 <- location.est
		location.est <- matrix(ncol=2,nrow=nn)
		for(i in 1:nn) {
			location.est[i,] <- c(location2[schlussel[1,i]],location2[schlussel[2,i]])
			}
		}

	# calculation of pairwise SSCM

	rho <- matrix(nrow=nn,ncol=2)
	for (i in 1:nn) {
		X2 <- cbind(X[,schlussel[1,i]],X[,schlussel[2,i]])
		X2 <- t(t(X2)-location.est[i,])
		Norm <- sqrt(X2[,1]^2+X2[,2]^2)
		index <- (Norm>0)
		Norm <- Norm[index]
		X2 <- X2[index,]
		ni <- sum(index)
		X2[,1] <- X2[,1]/Norm
		X2[,2] <- X2[,2]/Norm
		rho[i,] <- (t(X2)%*%X2/ni)[1,]
	}

	# calculation of correlation coefficients
	rho <- transformationcor(rho[,1],rho[,2]) 
	erg <- matrix(0,ncol=p,nrow=p)

	for (i in 1:(p-1)) {
		for (j in (i+1):p)	{
				erg[i,j] <- rho[(i-1)*p-(i-1)*i/2+j-i]
			}
		}
	erg <- erg+t(erg)
	diag(erg) <- 1

	return(erg)
	}

X <- t(t(X) - location.est)
Norm <- sqrt(apply(X^2,1,sum))
index <- Norm>0
ni <- sum(index)
X <- X[index,]
Norm <- Norm[index]
X <- X/Norm
sscm <- (t(X)%*%X/ni)
return(cov2cor(SSCM2Shape(sscm,...)))
}



