`fthetaf` <-
function(angle,		# Leaf angle
				 angledistobj=NULL,
				 degrees=FALSE,	# If TRUE, leaf angle is given in degrees, otherwise in radians.
				 distribution,
				 distpars=NA
 
				 ){
    
	if(!is.null(angledistobj)){
		distribution <- angledistobj$distribution
		distpars <- angledistobj$distpars
	}
	
	if(degrees)angle <- angle * pi/180
	
	if(distribution %in% c("ellipsoid","ellipsoidal")){
	
		if(all(is.na(distpars)))stop("Must provide X parameter\n")
		X <- distpars[1]
	
		# Spherical leaf angle distribution:
		if(X == 1){
			res <- sin(angle)
		} else {
		
    		# Wang et al. 2007 (AgForMet)
    		if(X < 1){
    			eps <- sqrt(1-X^2)
    			lambda <- X + asin(eps) / eps
    		}
    		if(X > 1){
    			eps <- sqrt(1 - X^-2)
    			lambda <- X + log((1+eps)/(1-eps))/(2*eps*X)
    		}
    		# Approximation: lambda <- X + 1.744*(X + 1.182)^-0.733
    		res <- 2 * X^3 * sin(angle) / (lambda * ( cos(angle)^2 + X^2*sin(angle)^2 )^2)	
        }
        
	}
	
	if(distribution == "rotatedell"){
	
		if(all(is.na(distpars)))stop("Must provide X parameter\n")
		X <- distpars[1]
	
		# Spherical leaf angle distribution:
		if(X == 1){
			res <- sin(angle)
		}
		if(X < 1){
			eps <- sqrt(1-X^2)
			lambda <- X + asin(eps) / eps
		}
		if(X > 1){
			eps <- sqrt(1 - X^-2)
			lambda <- X + log((1+eps)/(1-eps))/(2*eps*X)
		}
		res <- 2 * X^3 * cos(angle) / (lambda * ( sin(angle)^2 + X^2*cos(angle)^2 )^2)
	}
	
	# de Wit's 6 functions.
	if(distribution == "planophile") res <- (2 / pi) * (1 - cos(2*angle))
	if(distribution == "erectophile") res <- (2 / pi) * (1 + cos(2*angle))
	if(distribution == "plagiophile") res <- (2 / pi) * (1 - cos(4*angle))
	if(distribution == "extremophile") res <- (2 / pi) * (1 + cos(4*angle))
	if(distribution == "spherical") res <- sin(angle)
	if(distribution == "uniform") res <- 2/pi

	# # Two-parameter Beta
	if(distribution == "twoparbeta"){
	
		if(all(is.na(distpars)))stop("Must provide distpars\n")
		alphamean <- distpars[1]
		tvar <- distpars[2]
	
		# Find parameters m and v of the beta distribution based on mean and variance of leaf angle:
		# Wang et al. 2007, Eqs. 23-26.
	
		tmean <- 2 * alphamean * (pi/180) / pi
		sigma0 <- tmean * (1-tmean)
		v <- tmean * (sigma0/tvar - 1)
		m <- (1-tmean) * (sigma0/tvar - 1)
		
		te <- 2 * angle / pi
		te <- ifelse(te == 0, 1E-09, te)
		te <- ifelse(te == 1, 1 - 1E-09, te)
		
		res <- dbeta(te, v, m) / (pi/2)
	}
	
	if(degrees)
		return(res * pi/180)
	else
		return(res)
}

