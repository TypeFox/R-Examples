dN.f <-
function(x, mw, sd, ...){#		pdf(N(mw,sd)) = 1/(sqrt(2*pi*sd2))*exp(-1/2*((x-mw)/sd)2)
#   ------------------------------------------------------------------------------------------------
			return(as.numeric(dnorm(x = x, mean = mw, sd = sd, log = FALSE)))
	}
