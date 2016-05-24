srm_gif <-
function(data, evalpts, params, TT = NA, tplus = FALSE)
{
	#  with the same algorithm as linksrm_gif
	if (is.null(data)) data <- cbind(time = Inf, magnitude = 0)
	times <- data[, "time"]
	magnitude <- data[, "magnitude"]
	A <- params[1]
	B <- params[2]
	cc <- params[3]
	if (tplus) FUN <- "<="
		else FUN <- "<"
	if(any(is.na(TT))) {
		if(is.vector(evalpts))
			eval.times <- evalpts
		else eval.times <- evalpts[, "time"]
		St <- matrix(10^(0.75 * magnitude), nrow = 1) %*% 
				outer(times, eval.times, FUN = FUN)	
		ci <- exp(A + B * (eval.times - cc*St))
	}
	else {
		within <- (times < TT[2]) & (times >= TT[1])
		ti <- c(times[within], TT[2])
		nt <- length(ti)
		St <- matrix(10^(0.75 * magnitude), nrow = 1) %*% 
			outer(times, ti, FUN = FUN)
		ci <- sum(exp(A - B*cc*St)/B * (exp(B*ti) - 
			exp(B*c(TT[1], ti[ - nt]))))
	}
	names(ci) <- NULL
	ci <- as.vector(ci)
	return(ci)
}
attr(srm_gif, "rate") <- "increasing"


