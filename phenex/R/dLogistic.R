.dLogistic <- 
function(ndvi){
	days <- length(ndvi)
	ndvi <- .linIP(ndvi)
	time <- 1:days

	Wx <- function(x){
		erg <- sum(((x[1] + (x[3]/(1+exp(-x[6]*(time-x[4])))) - 
			((x[3]+x[1]-x[2])/(1+exp(-x[5]*(time-x[7]))))) - (ndvi))^2)
		return(ifelse((is.infinite(erg)||is.nan(erg)),1e50,erg))
	}

	Wt <- function(vb,ve,k,p,d,c,q) {
   		erg <- vb + (k/(1+exp(-c*(time-p)))) - ((k+vb-ve)/(1+exp(-d*(time-q))))
		return(erg)
  	}
	
	vb <- ve <- 0.1
    	c <- d <- 0.1
    	k <- 0.6
    	p <- 90
	q <- 320

	optimal <- optim(c(vb,ve,k,p,d,c,q),fn=Wx,gr=NULL,method="L-BFGS-B",
			lower=c(0.001,0.001,0.001,1,0.001,0.001,1), upper=c(0.3,0.3,1,300,0.5,0.5,340),
			control=list(maxit = 5000, pgtol = 1e-10, 
			ndeps = c(1e-10, 1e-10, 1e-10, 1, 1e-10, 1e-10,1), 
			lmm = 200))	

	vb <- optimal$par[1]
	ve <- optimal$par[2]
  	k <- optimal$par[3]
	p <- optimal$par[4]
  	d <- optimal$par[5]
	c <- optimal$par[6]
  	q <- optimal$par[7]

	model.optim <- Wt(vb,ve,k,p,d,c,q)

	count <- 10
	tolerance <- 1e-5
	repeat {
		model <- try(nls(ndvi ~ Wt(vb,ve,k,p,d,c,q),
                		start=list(vb=vb,ve=ve,k=k,p=p,d=d,c=c,q=q), 
				control=list(maxiter=200,tol=tolerance,
				minFactor=1/4096)), silent=TRUE)
		if (inherits(model, "try-error")==FALSE){
			model.interpol <- predict(model)
			break
		} else {
			count <- count - 1
			if (count==0){
				return(rep(NA, days))
			}
			tolerance <- tolerance*10
		}
	}

	return(model.interpol)
}
