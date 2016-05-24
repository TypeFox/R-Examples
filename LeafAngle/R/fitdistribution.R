`fitdistribution` <-
function(angles, 				# vector of angles (radians or degrees, see degrees argument).
distribution, 			# distribution to be fit (character string)
fitmethod=c("loglik","chisq"),
ellipsmethod=2, 		# method for finding X in the ellipsoidal. 1= Wang method, 2=optimize chisq or loglik.
degrees=TRUE,			# if FALSE, input 'angles' is in radians, otherwise in degrees.
...){

	fitmethod <- match.arg(fitmethod)
							
	# Convert degrees properly:
	if(degrees){
		anglesrad <- angles * pi/180
	}
	if(!degrees){
		anglesrad <- angles
		angles <- anglesrad * 180/pi
	}

    # Replace 0 with a small number, and 90 with 90 - a small number.
    angles[angles < 0.01] <- 0.01
    angles[angles > (90-0.01)] <- 90 - 0.01

	# Toss out angles > 90, print warning.
	nover <- length(angles[angles > 90])
	
	if(nover > 0){
		warning("Deleting ",nover," angles > 90degrees")
		angles <- angles[angles <=90]
		anglesrad <- anglesrad[angles <=90]
    }
	
	# Ellipsoidal and rotated ellipsoidal
	if(distribution %in% c("ellipsoid","rotatedell")){
	
		if(fitmethod == "chisq"){
			# Eq. 30
			if(ellipsmethod == 1)Xfit <- -3 + (mean(anglesrad) / 9.65)^-0.6061
			if(ellipsmethod == 2){
				chisq_wrapper <- function(Xs, angles)chisqfit(angles, distribution=distribution, distpars=Xs)
				Xfit <- optimize(chisq_wrapper, interval=c(0.05,20), angles)$minimum
			}
			chisq <- chisqfit(angles, distribution=distribution, distpars=Xfit)
			
			l <- list(distribution=distribution,chisq=chisq, loglik=NA, AIC=NA,  distpars=Xfit,
						fitmethod="Chi-squared", fitdata = angles)
			class(l) <- "angledist"
			return(l)
		}
		if(fitmethod == "loglik"){
			loglik_wrapper <- function(Xs, angles)-loglikfit(angles, distribution=distribution, distpars=Xs)
			Xfit <- optimize(loglik_wrapper, interval=c(0.05,20), angles)$minimum
			loglik <- loglikfit(angles, distribution=distribution, distpars=Xfit)
			AICfit <- -2*loglik + 2*1  # one parameter
			l <- list(distribution=distribution,chisq=NA, loglik=loglik, AIC=AICfit, distpars=Xfit,
						fitmethod="Log-likelihood", fitdata = angles)
			class(l) <- "angledist"
			return(l)
		}
	}

	# Beta (two-parameter)
	if(distribution == "twoparbeta"){
		
		tvals <- 2 * anglesrad / pi
		tvarfit <- var(tvals)
		alphameanfit <- mean(angles)
	
		if(fitmethod == "chisq"){
			chisq <- chisqfit(angles, distribution="twoparbeta", distpars=c(alphameanfit,tvarfit))
			
			l <- list(distribution="twoparbeta", chisq=chisq, loglik=NA, AIC=NA, distpars=c(alphameanfit,tvarfit),
						fitmethod="Chi-squared", fitdata = angles)
			class(l) <- "angledist"
			return(l)
		}
		if(fitmethod == "loglik"){
			loglik <- loglikfit(angles, distribution="twoparbeta", distpars=c(alphameanfit,tvarfit))
			AICfit <- -2*loglik + 2*2  # two parameters
			l <- list(distribution="twoparbeta", chisq=NA, loglik=loglik, AIC=AICfit, distpars=c(alphameanfit,tvarfit),
						fitmethod="Log-likelihood", fitdata = angles)
			class(l) <- "angledist"
			return(l)
		}
	}
	
	# One of de Wit's functions (no parameters):
	if(distribution %in% c("planophile","erectophile","plagiophile","extremophile","spherical","uniform")){
		
		if(fitmethod=="chisq"){
			chisq <- chisqfit(angles, distribution=distribution)
			l <- list(distribution=distribution,chisq=chisq, loglik=NA, AIC=NA, distpars=NA,
						fitmethod="Chi-squared", fitdata = angles)
			class(l) <- "angledist"
			return(l)
		}
		if(fitmethod=="loglik"){
			loglik <- loglikfit(angles, distribution=distribution)
			AICfit <- -2*loglik + 2*0  # no parameters
			l <- list(distribution=distribution,chisq=NA, loglik=loglik, AIC=AICfit, distpars=NA,
						fitmethod="Log-likelihood", fitdata = angles)
			class(l) <- "angledist"
			return(l)
		}
	}

# Should have returned already:
stop("Unknown distribution")
}

