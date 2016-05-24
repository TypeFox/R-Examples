f.grovest3 <- function(formula, data, mu = T, tau = T, covar.est = F, init.plot = T)
{
#
# COMPUTE STARTING VALUES USING EXPLICIT MLE FOR SOME GROUPS/SUBGROUPS
#
#
#
	f.est <- function(response, mu, tau, init.plot = T)
	{
#
# FIND MLE FOR mu AND c (DISREGARDING COVARIATES), COMPUTED FROM THE KAPLAN-MEIER
# ESTIMATE OF THE SURVIVAL FUNCTION. IF mu IS GIVEN A NUMERIC VALUE ONLY THE VALUE OF 
# c IS COMPUTED WHILST mu IS KEPT FIXED.
#
		fit <- survfit(formula = response ~ 1)	#

		surv <- fit$surv
		time <- fit$time	#
#
# ESTIMATE MEAN AND HARMONIC MEAN
		est.mean <- time[1] + sum(diff(time) * surv[ - length(surv)])	#
		est.mean.1dT <- sum((1 - surv[ - length(surv)]) * ( - diff(1/
			time)))
		est.harm.mean <- 1/est.mean.1dT	#
#
		if(is.numeric(mu)) {
# ESTIMATE c, KEEP mu FIXED
			est.mu <- mu
			est.c <- (mu * est.harm.mean + sqrt(mu^2 * 
				est.harm.mean^2 + 4 * est.harm.mean))/2
		}
		else {
# ESTIMATE BOTH mu AND c
			est.mu <- sqrt(est.harm.mean/(est.mean^2 - 
				est.harm.mean * est.mean))
			est.c <- est.mu * est.mean
		}
		if(is.numeric(tau))
			est.tau <- tau
		else est.tau <- est.mu
		if(init.plot) {
			plot(fit, conf.int = F, main = 
				"Compare KM and initial estimate of param. surv."
				)
			lines(time, f.B(time, mu = est.mu, tau = est.tau, c1 = 
				est.c))
		}
#
		return(c(est.mu = est.mu, est.tau = est.tau, est.c = est.c))
	}
#
# EXTRACT time, status AND covariates
model <- model.frame(formula, data = data)
response <- model.extract(model, "response")
est <- f.est(response = response, mu = mu, tau = tau, init.plot = 
	init.plot)
cat("Initial estimates:\n")
print(est)	#
cat("*************\n")	#
return(est)	#
#
#
# FUNCTION ENDS HERE, THE REST IS INACTIVE (FOR THE TIME BEING)
#
#
## Folgende kode er kommentert vekk for aa slippe feilmeld. i CRAN,
## men kan evt. aktiveres senere
##
#	if(covar.est) {
#		deleverdier <- c(sex = 1, cond = 1, tstage = 2, nstage = 1)
#		for(kovartext in c("sex", "cond", "tstage", "nstage")) {
#			#
#			kovar <- data[[kovartext]]	#
#			under.deling <- kovar <= deleverdier[kovartext]	#
#			skaler <- mean(kovar[!under.deling]) - mean(kovar[
#				under.deling])
#			kovar[under.deling] <- 1
#			kovar[!under.deling] <- 2
#			assign("kovar", kovar, f = 1)
#			par(mfcol = c(2, 1))
#			cat("covariatestimering for ", kovartext, "\n")
#			est.lower <- f.est(response = response[under.deling,  ],
#				mu = est.mu)
#			est.upper <- f.est(response = response[!under.deling,  
#				], mu = est.mu)	#
#			#  plot(fit, conf.int = F, main = kovartext)
#			#  lines(fit$time, f.B(t = fit$time, mu = mu, c1 = est.c[i]), lty # = 2) #
#			#  cat(kovartext, "\n")#
#			#   f.vis((est.upper["est.c"] - est.lower["est.c"])/skaler)
#		}
#	}
}

