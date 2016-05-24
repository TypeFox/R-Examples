.dSig <-
function(ndvi){
	days <- length(ndvi)
	ndvi <- .linIP(ndvi)
	time <- 1:days

	DSx <- function(x){
		pos1 <- x[1]
		pos2 <- x[2]
		width1 <- x[3]
		width2 <- x[4]
		erg <- sum(((0.5*(tanh((time-pos1)/width1)-
				tanh((time-pos2)/width2)))-ndvi)^2)
		return(ifelse((is.infinite(erg)||is.nan(erg)),1e50,erg))
	}

	maxind <- order(ndvi, decreasing=TRUE)[1]
	#calculate double sigmoid function
	pos1 <- 0.75*maxind
	pos2 <- 1.25*maxind
	width1 <- maxind
	width2 <- days-maxind

	optimal <- optim(c(pos1,pos2,width1,width2),fn=DSx,gr=NULL,method="L-BFGS-B",
				lower=c(0,0,0,0), upper=c(365,365,365,365),
				control=list(maxit = 5000, pgtol = 1e-10, 
				ndeps = c(1, 1, 1, 1), lmm = 200))
	pos1 <- optimal$par[1]
	pos2 <- optimal$par[2]
	width1 <- optimal$par[3]
	width2 <- optimal$par[4]  

	count <- 10
	tolerance <- 1e-5
	repeat {
		model.nls <- try(nls(ndvi ~ 
			0.5*(tanh((time-pos1)/width1)-tanh((time-pos2)/width2)), 
			start=list(pos1=pos1, pos2=pos2, width1=width1, width2=width2),
			control=list(maxiter=200, tol=tolerance,minFactor=1/4096)), 
			silent=TRUE)
		if (inherits(model.nls, "try-error")==FALSE){
			model <- predict(model.nls)
			break
		} else {
			count <- count-1
			if (count==0){
				return(rep(NA, days))
			}
			tolerance <- tolerance*10
		}
	}

	return(model)
}
