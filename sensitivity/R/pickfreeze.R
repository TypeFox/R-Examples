sobolpickfreeze <- function(y1,y2,nboot){
	output <- (apply(y1*y2,1,mean) - apply(y1,1,mean)*apply(y2,1,mean))/apply(y1,1,var)
	if(nboot > 1){
		n <- dim(y1)[2]
		for (i in 1:nboot){
			b <- sample(n,replace = TRUE)
			output <- rbind(output,(apply(y1[,b]*y2[,b],1,mean) - apply(y1[,b],1,mean)*apply(y2[,b],1,mean))/apply(y1[,b],1,var))
		}
	}
	return(output)
}

soboljansenpickfreeze <- function(y1,y2,nboot){
	varY <- apply(cbind(y1,y2),1,var)
	output <- (varY-apply((y1-y2)^2,1,mean)/2)/varY
	if(nboot > 1){
		n <- dim(y1)[2]
		for (i in 1:nboot){
			b <- sample(n,replace = TRUE)
			varY <- apply(cbind(y1[,b],y2[,b]),1,var)		
			output <- rbind(output,(varY-apply((y1[,b]-y2[,b])^2,1,mean)/2)/varY)
		}
	}
	return(output)
}

sobolEffpickfreeze <- function(y1,y2,nboot){
	output <- (apply(y1*y2,1,mean) - apply(cbind(y1,y2),1,mean)^2)/apply(cbind(y1,y2),1,var)
	if(nboot > 1){
		n <- dim(y1)[2]
		for (i in 1:nboot){
			b <- sample(n,replace = TRUE)
			output <- rbind(output,(apply(y1[,b]*y2[,b],1,mean) - apply(cbind(y1[,b],y2[,b]),1,mean)^2)/apply(cbind(y1[,b],y2[,b]),1,var))
		}
	}
	return(output)
}

sobolT2002pickfreeze <- function(y1,y2,nboot){
	output <- (1-(apply(y1*y2,1,mean) - apply(y1,1,mean)^2)/apply(y1,1,var))
	if(nboot > 1){
		n <- dim(y1)[2]
		for (i in 1:nboot){
			b <- sample(n,replace = TRUE)
			output <- rbind(output,(1-(apply(y1[,b]*y2[,b],1,mean) - apply(y1[,b],1,mean)^2)/apply(y1[,b],1,var)))
		}
	}
	return(output)
}

sobolTjansenpickfreeze <- function(y1,y2,nboot){
	output <- (apply((y1-y2)^2,1,mean)/2)/apply(y1,1,var)
	if(nboot > 1){
		n <- dim(y1)[2]
		for (i in 1:nboot){
			b <- sample(n,replace = TRUE)
			output <- rbind(output,(apply((y1[,b]-y2[,b])^2,1,mean)/2)/apply(y1[,b],1,var))
		}
	}
	return(output)
}

sobol2007pickfreeze <- function(y1,y2,y3,nboot){
	output <- apply(y3*(y2-y1),1,mean)/apply(y1,1,var)
	if(nboot > 1){
		n <- dim(y1)[2]
		for (i in 1:nboot){
			b <- sample(n,replace = TRUE)
			output <- rbind(output,apply(y3[,b]*(y2[,b]-y1[,b]),1,mean)/apply(y1[,b],1,var))
		}
	}
	return(output)
}

sobolT2007pickfreeze <- function(y1,y2,nboot){
	output <- apply(y1*(y1-y2),1,mean)/apply(y1,1,var)
	if(nboot > 1){
		n <- dim(y1)[2]
		for (i in 1:nboot){
			b <- sample(n,replace = TRUE)
			output <- rbind(output,apply(y1[,b]*(y1[,b]-y2[,b]),1,mean)/apply(y1[,b],1,var))
		}
	}
	return(output)
}












