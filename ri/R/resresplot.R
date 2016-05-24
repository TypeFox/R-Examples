resresplot <-
function(Y,Z,X,prob=NULL,scale=1) {
		if(is.null(prob)) {
		warning("Probabilities not specified. Assuming equal probabilities.")
		prob <- rep(.5,length(Y))
	}
	w = Z/prob + (1-Z)/(1-prob)
	w = w/mean(w)
	Yres <- Y - lm(Y~X,weights=w)$fitted
	Zres <- Z - lm(Z~X,weights=w)$fitted
	reg <- lm(Yres~Zres,weights=w)
	symbols(Zres,Yres,circles=scale*.1*w/diff(range(Yres)),inches=FALSE,xlab="Residualized Z",ylab="Residualized Y",main=paste("Estimated ATE = ",round(reg$coefficients["Zres"],3)))
	abline(reg)
	}
