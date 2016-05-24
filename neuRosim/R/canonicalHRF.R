canonicalHRF <-
function(x, param=NULL, verbose=TRUE){
	if(is.null(param)){
		if(verbose==TRUE){
			warning("Default parameters for HRF are used")
		}
        	param <- list()
		param$a1 <- 6
		param$a2 <- 12
		param$b1 <- 0.9
		param$b2 <- 0.9
		param$c <- 0.35
	}
        d1 <- param$a1*param$b1		# time to peak of response
        d2 <- param$a2*param$b2		# time to peak of undershoot

        (x/d1)^param$a1*exp(-(x-d1)/param$b1) - param$c*(x/d2)^param$a2*exp(-(x-d2)/param$b2)
}

