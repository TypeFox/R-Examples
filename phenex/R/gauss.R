.gauss <-
function(ndvi, asym=FALSE){
	if (!is.logical(asym)){ stop("'asym' should be of type 'logical'") }
	days <- length(ndvi)
	model <- vector(mode="numeric",length=days)
	ndvi.interpol <- .linIP(ndvi)

	#parameters of gaussian function
	maxpos <- order(ndvi, decreasing=TRUE)
	maxpos <- maxpos[which(!is.na(ndvi[maxpos]))]
	chosen <- which(maxpos <= 210 & maxpos >= 120)
	if (length(chosen)==0){
		mu <- maxpos[1]
	} else {
		mu <- maxpos[chosen[1]]
	}
	
	sig <- sqrt(var(na.omit(ndvi.interpol)))
	base <- mean(ndvi.interpol[1:(days/6)])
	scal <- 1
	time <- 1:days

	#calculate gaussian function
	if (asym){
		ndvi <- ifelse(is.na(ndvi), -1, ndvi)
		model <- .C("asymgauss", rdays=as.integer(days), ndvi=as.numeric(ndvi),
			mustart=as.integer(mu), sigstart=as.numeric(sig), 
			rbase=as.numeric(base), model=as.numeric(model), 
			PACKAGE="phenex")$model
	} else {
		Gx <- function(x){
			mu <- x[1]
			sig <- x[2]
			scal <- x[3]
			erg <- sum(((((scal / (sig * sqrt(2 * pi)))*exp(-0.5*((((time/days) - 
					(mu/days)) / sig)^2)))+base)-ndvi.interpol)^2)
			return(ifelse((is.infinite(erg)||is.nan(erg)),1e50,erg))
		}

		G <- function(mu, sig, scal, base, time, days){
			erg <- ((scal / (sig * sqrt(2 * pi)))*
					exp(-0.5*((((time/rep(days,days)) - 
					(mu/days)) / sig)^2)))+base
			return(erg)
		} 

		optimal <- optim(c(mu,sig,scal),fn=Gx,gr=NULL,method="L-BFGS-B",
					lower=c(120.0,0.001,0.1), upper=c(210.0,2.0,2.0),
					control=list(maxit = 5000, pgtol = 1e-10, 
					ndeps = c(1e-10, 1e-10, 1e-10), lmm = 200))
		mu <- optimal$par[1]
		sig <- optimal$par[2]
		scal <- optimal$par[3]

		count <- 10
		tolerance <- 1e-5
		repeat {
			model.nls <- try(nls(ndvi.interpol ~ G(mu=mu, sig=sig, scal=scal, 
						base=base, time=time, days=days), 
						start=list(sig=sig, scal=scal, mu=mu), 
						control=list(maxiter=200, tol=tolerance,
						minFactor=1/4096)),silent=TRUE)
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
	}

	return(model)
}
