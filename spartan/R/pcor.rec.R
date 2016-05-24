pcor.rec <-
function(x,y,z,corMethod="p",na.rm=TRUE){

	x <- c(x)
	y <- c(y)
	z <- as.data.frame(z,check.names=FALSE)

	if(dim(z)[2] == 0){
		stop("There should be given data\n")
	}

	data <- data.frame(x,y,z,check.names=FALSE)

	if(na.rm == TRUE){
		data = na.omit(data)
	}

	# recursive formula
	if(dim(z)[2] == 1){
		tdata <- na.omit(data.frame(data[,1],data[,2]),check.names=FALSE)
		rxy <- cor(tdata[,1],tdata[,2],corMethod)

		tdata <- na.omit(data.frame(data[,1],data[,-c(1,2)]),check.names=FALSE)
		rxz <- cor(tdata[,1],tdata[,2],corMethod)

		tdata <- na.omit(data.frame(data[,2],data[,-c(1,2)]),check.names=FALSE)
		ryz <- cor(tdata[,1],tdata[,2],corMethod)

		rxy.z <- (rxy - rxz*ryz)/( sqrt(1-rxz^2)*sqrt(1-ryz^2) )
		
		return(rxy.z)
	}else{
		x <- c(data[,1])
		y <- c(data[,2])
		z0 <- c(data[,3])
		zc <- as.data.frame(data[,-c(1,2,3)],check.names=FALSE)

		rxy.zc <- pcor.rec(x,y,zc,corMethod=corMethod,na.rm=na.rm)
		rxz0.zc <- pcor.rec(x,z0,zc,corMethod=corMethod,na.rm=na.rm)
		ryz0.zc <- pcor.rec(y,z0,zc,corMethod=corMethod,na.rm=na.rm)
		
		rxy.z <- (rxy.zc - rxz0.zc*ryz0.zc)/( sqrt(1-rxz0.zc^2)*sqrt(1-ryz0.zc^2) )
		return(rxy.z)
	}			
}

