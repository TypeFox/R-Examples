reco <- function(R, Temp, Tref = 10, T0 = -46.02, method = "all", min.dp = 6){
	## transform temperatures
#	Temp <- Temp + 273.15
#	Tref <- Tref + 273.15
#	T0 <- T0 + 273.15
	## check number of cm
	if(length(R) < min.dp){
		stop("Not enough concentration measurements")
	}
	## linear model
	res <- list(linear = lm(R ~ Temp))
	## prepare start value for Rref (global, t1 here)
	t1 <- coef(res$linear)[1]
	## non-linear models
	# Arr
	tryArr <- function(E0){
		try(nls(R ~ t1 * exp(E0 * (1/(Tref-T0) - 1/(Temp-T0))), start=list(t1 = t1, E0 = E0), model = TRUE), silent = TRUE)
		}
	Arr <- lapply(exp(seq(0,10,0.5)), function(x) tryArr(x))
	if(sum(sapply(Arr, class) != "try-error")==0){
		Arr <- NA
		class(Arr) <- "try-error"
		t2 = 10
	}
	else{
		Arr <- Arr[sapply(Arr, class) != "try-error"]
		Arr <- Arr[which.min(sapply(Arr, function(x) x[]$convInfo$finIter))][[1]]
		t2 <- coef(Arr)[2]
	}
	t3 <- t2/2
	# Q10
	Q10 <- try(nls(R ~ t1 * t2^((Temp-Tref)/10), start=list(t1 = 1, t2 = 1), model = TRUE), silent = TRUE)
	# LT
	LT <- try(nls(R ~ t1 * exp(-t2 / (Temp+273.15-t3)), start=list(t1 = t1, t2 = 308.56, t3 = 227.13), model = TRUE, control=nls.control(minFactor=1/4096)), silent = TRUE)
	# restricted LT
	LTR <- try(nls(R ~ t1 * exp(-308.56 / (Temp+46.02)), start=list(t1 = t1), model = TRUE), silent = TRUE)
	# Logistic
	Log <- try(nls(R ~ t1 / (1 + exp(t2-t3*Temp)), start=list(t1 = t1, t2 = t2, t3 = t3), model = TRUE), silent = TRUE)
	
	## compile results
	res$Arr <- Arr
	res$Q10 <- Q10
	res$LT <- LT
	res$LTR <- LTR
	res$Log <- Log
	
	## prepare return
	METHODS <- c("all", "not.failed", "best", "linear", "arrhenius", "Q10", "lt", "lt.restricted", "logistic")
	method <- pmatch(method, METHODS)
	
	if(method == 2){
		## sort out the failed attempts
		nerrs <- lapply(res, class) != "try-error"
		res <- res[nerrs]
	}
	if(method == 3){
		## sort out the failed attempts
		nerrs <- lapply(res, class) != "try-error"
		res <- res[nerrs]		
		res <- res[which.min(sapply(res, AIC))]
	}
	if(method == 4){
		res <- res[1]
	}
	if(method == 5){
		res <- res[2]
	}
	if(method == 6){
		res <- res[3]
	}
	if(method == 7){
		res <- res[4]
	}
	if(method == 8){
		res <- res[5]
	}	
	if(method == 9){
		res <- res[6]
	}					
	class(res) <- "reco"
	return(res)
}
