predict.MuFicokm <- function (object, newdata, type, se.compute=TRUE, 
                                  cov.compute=FALSE, checkNames=FALSE, ...)
{
	nlevel <- object$nlevel
	## initialisation
	if(as.numeric(length(type))==1){
		type1 <- type
	}else{
		type1 <- type[[1]]
	}
	if(as.numeric(length(se.compute))==1){
		se.compute1 <- se.compute
	}else{
		se.compute1 <- se.compute[[1]]
	}
	if(as.numeric(length(cov.compute))==1){
		cov.compute1 <- cov.compute
	}else{
		cov.compute1 <- cov.compute[[1]]
	}
	if(as.numeric(length(checkNames))==1){
		checkNames1 <- checkNames
	}else{
		checkNames1 <- checkNames[[1]]
	}

	object1 <- object$cok[[1]]
	z1x <- predict(object = object1,
				newdata = data.frame(newdata),
				type = type1, 
				se.compute=se.compute1, 
                        cov.compute=cov.compute1, 
				checkNames=checkNames1)
	mux <- list()
	mux[[1]] <- z1x$mean
	varx <- list()
	varx[[1]] <-  z1x$sd^2
	if(identical(as.numeric(length(cov.compute)),1)){
		if(identical(cov.compute,TRUE)){
			CovMat <- list()
			CovMat[[1]] <- z1x$cov
		}
	}

	for(i in 2:nlevel){

		## initialisation
		if(as.numeric(length(type))==1){
			typei <- type
		}else{
			typei <- type[[i]]
		}
		if(as.numeric(length(se.compute))==1){
			se.computei <- se.compute
		}else{
			se.computei <- se.compute[[i]]
		}
		if(as.numeric(length(cov.compute))==1){
			cov.computei <- cov.compute
		}else{
			cov.computei <- cov.compute[[i]]
		}
		if(as.numeric(length(checkNames))==1){
			checkNamesi <- checkNames
		}else{
			checkNamesi <- checkNames[[i]]
		}

		objecti <- object$cok[[i]]
		zix <- predict.kmCok(
					object = objecti,
					newdata = data.frame(newdata),
					newZ = mux[[i-1]],
					type = typei, 
					se.compute=se.computei, 
                              cov.compute=cov.computei, 
					checkNames=checkNamesi)
		mux[[i]] <- zix$mean
		varx[[i]] <- objecti@trend.coef[1]^2*varx[[i-1]] + zix$sd^2

		if(identical(as.numeric(length(cov.compute)),1)){
			if(identical(cov.compute,TRUE)){
				CovMat[[i]] <- objecti@trend.coef[1]^2*CovMat[[i-1]] + zix$cov
			}
		}
		
	}
	if(identical(as.numeric(length(cov.compute)),1)){
		if(identical(cov.compute,TRUE)){
			output <- list(mean=mux[[nlevel]],sig2=varx[[nlevel]],C=CovMat[[nlevel]],mux = mux,varx=varx,CovMat=CovMat)
		}else{
			output <- list(mean=mux[[nlevel]],sig2=varx[[nlevel]],mux = mux,varx=varx)
		}
	}else{
		output <- list(mean=mux[[nlevel]],sig2=varx[[nlevel]],mux = mux,varx=varx)
	}

	return(output)
}

