.growth <- 
function(ndvi){
	days <- length(ndvi)
	ndvi <- .linIP(ndvi)
	
	time <- 1:days
	Wx <- function(x){
		W <- x[1]
		a <- x[2] 
		r <- x[3]
		p <- x[4] 
		my <- x[5]
		p <- log(abs(p))
		erg <- sum(((W * ((a+1)^(r/exp(p))) * exp(((-1) * my) * time)) / 
				((1+a*exp(-exp(p)*time))^(r/exp(p)))-ndvi)^2)
		return(ifelse((is.infinite(erg)||is.nan(erg)),1e50,erg))
	}

	Wt <- function(W, a, r, p, my, time){
		p <- log(abs(p))
		erg <- (W * ((a+1)^(r/exp(p))) * exp(((-1) * my) * time)) / 
			((1+a*exp(-exp(p)*time))^(r/exp(p)))
		return(erg)
	}
	
	W.base <- mean(na.omit(ndvi[1:(days/5)]))

	tmax <- order(ndvi, decreasing=TRUE)[1]
	r <- ((((max(na.omit(ndvi))/W.base)^(1 / as.numeric(tmax)))) - 1)*6
	p <- r
	my <- r/3
	a <- (exp(tmax*p) / ((r/my)-1))
	time <- 1:days

	optimal <- optim(c(W.base,a,r,p,my),fn=Wx,gr=NULL,method="L-BFGS-B",
				lower=c(0.0,0.001,0.001,0.001,0.005), 
				upper=c(0.3,1e4,0.05,3.0,0.05),
				control=list(maxit = 5000, pgtol = 1e-10, 
				ndeps = c(1e-10, 1e-10, 1e-10, 1e-10, 1e-10), 
				lmm = 200))
	W.base <- optimal$par[1]
	a <- optimal$par[2]
	r <- optimal$par[3]
	p <- optimal$par[4]
	my <- optimal$par[5]

	count <- 10
	tolerance <- 1e-5
	repeat {
		model.nls <- try(nls(ndvi ~ Wt(W,a,r,p,my,time), 
				start=list(W=W.base, a=a, r=r, p=p, my=my), 
				control=list(maxiter=200, tol=tolerance, minFactor=1/4096)),
				silent=TRUE)
		if (inherits(model.nls, "try-error")==FALSE){
			model <- predict(model.nls)
			break
		} else {
			count <- count - 1
			if (count==0){
				return(rep(NA, days))
			}
			tolerance <- tolerance*10
		}
	}

	return(model)
}
