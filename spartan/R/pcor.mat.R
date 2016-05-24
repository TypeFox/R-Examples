pcor.mat <-
function(x,y,z,corMethod="p",na.rm=TRUE){

	x <- c(x)
	y <- c(y)
	z <- as.data.frame(z,check.names=FALSE)

	if(dim(z)[2] == 0){
		stop("There should be given data\n")
	}

	data <- data.frame(x,y,z, check.names=FALSE)

	if(na.rm == TRUE){
		data = na.omit(data)
	}

	xdata <- na.omit(data.frame(data[,c(1,2)]),check.names=FALSE)
	Sxx <- cov(xdata,xdata,method=corMethod)

	xzdata <- na.omit(data)
	xdata <- data.frame(xzdata[,c(1,2)],check.names=FALSE)
	zdata <- data.frame(xzdata[,-c(1,2)],check.names=FALSE)
	Sxz <- cov(xdata,zdata,method=corMethod)

	zdata <- na.omit(data.frame(data[,-c(1,2)]),check.names=FALSE)
	Szz <- cov(zdata,zdata,method=corMethod)

	# is Szz positive definite?
	zz.ev <- eigen(Szz)$values
	if(min(zz.ev)[1]<0){
		stop("\'Szz\' is not positive definite!\n")
	}

	# partial correlation
	Sxx.z <- Sxx - Sxz %*% solve(Szz) %*% t(Sxz)
	
	rxx.z <- cov2cor(Sxx.z)[1,2]

	rxx.z
}

