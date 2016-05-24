etas_gif <-
function(data, evalpts, params, TT = NA, tplus = FALSE)
{
#----------------------------------------------------------------------
# calculate the conditional intensity at given times as a function
# of past event history up to each time.
#
# data           = converted catalogue data (time, magnitude matrix)
# eval.times     = time instances at which conditional intensity is to be 
#                  calculated
# params         = vector of parameters, with names attributes 
#                  described below
#  mu            = background rate
#  A             = proportionality constant
#  alpha         = sensitivity coefficient
#  CC            = time constant
#  P             = decay coefficient
#  TT		 = time to evaluate, if not NA, only integral term evaluated
#----------------------------------------------------------------------
	mu <- params[1]
	A <- params[2]
	alpha <- params[3]
	CC <- params[4]
	P <- params[5]	
	#------------------------------------------------------------
	if (is.null(data)) data <- cbind(time = Inf, magnitude = 0)
	magnitudes <- data[, "magnitude"]
	times <- data[, "time"]
	if(any(is.na(TT))) {
		#----------------------------------------------------	
		#  Calculate the cumulative sum, S(t), say. Can only
		#  use apply if calculating conditional intensity for 
		#  more than one time otherwise simply calculate directly
		#  if TT is not missing, then just evaluate the integral term
		#---------------------------------------------------
		if(is.vector(evalpts)) eval.times <- evalpts
		else eval.times <- evalpts[, "time"]
		if(length(eval.times) > 1) {
			St <- apply(t(eval.times), 2, function(z, magnitudes, 
				times, mu, A, P, CC, alpha, tplus)
			{
				if (!tplus) use <- times < z
					else use <- times <= z
				if(sum(use) > 0)
				  sum(exp(alpha * magnitudes[use]) * (1 + (z - 
				    times[use])/CC)^( - P))
				else 0
			}
			, magnitudes = magnitudes, times = times, mu = mu,
				A = A, P = P, CC = CC, alpha = alpha,
				tplus = tplus)
		}
		else {
			if (!tplus) use <- times < eval.times
				else use <- times <= eval.times
			if(sum(use) > 0)
				St <- sum(exp(alpha * magnitudes[use]) * (1 + (
				  eval.times - times[use])/CC)^( - P))
			else St <- 0
		}
		ci <- mu + A * St
	}
	else {
		#  if TT is not missing, then just evaluate the integral term
		S <- c(NA, NA)
		for(Ii in 1:2) {
			Z0 <- (times < TT[Ii])
			if(any(Z0)) {
				if(P != 1)
				  S[Ii] <- sum(exp(alpha * data[, "magnitude"][
				    Z0]) * (1 - (1 + (TT[Ii] - times[Z0])/CC)^( - 
				    P + 1)))/(P - 1)
				else S[Ii] <- sum(exp(alpha * data[, "magnitude"
				    ][Z0]) * log(1 + (TT[Ii] - times[Z0])/CC))
			}
			else S[Ii] <- 0
		}
		ci <- mu * (TT[2] - TT[1]) + CC * A * (S[2] - S[1])
	}
	names(ci) <- NULL
	return(ci)
}
attr(etas_gif, "rate") <- "decreasing"

