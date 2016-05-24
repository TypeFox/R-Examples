### R code from vignette source 'pdmod.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: trials
###################################################
library(pdmod)
rewards = c(rep(1, 20), rep(0, 10))
schedule = c(1:10, (2*TV_DAY):(2*TV_DAY + 9), (4*TV_DAY):(4*TV_DAY + 1), 
			 (6*TV_DAY):(6*TV_DAY + 1), (8*TV_DAY):(8*TV_DAY + 1), 
			 (10*TV_DAY):(10*TV_DAY + 1), (12*TV_DAY):(12*TV_DAY + 1))
trials = TimedVector(rewards, schedule)
trials


###################################################
### code chunk number 2: compute_model
###################################################
params = c(0.9, 0.01, 0.04, 0.4, 0.25, 4.5, 500)
est = computeModel(trials, mFast=params[1], mSlow=params[2], 
				   n=params[3], h=params[4], g=params[7], 
				   verbose = TRUE)
plot(est, actual = rep(NA, length(trials)))


###################################################
### code chunk number 3: fit_model
###################################################
responses = c(0.001, 0.4, 1.4, 2.5, 3.6, 3.7, 3.6, 4.2, 4.6, 4.5, 
			  4.8, 5.0, 4.9, 4.8, 5.1, 4.9, 4.9, 5.0, 4.8, 4.9, 
			  3.7, 1.1, 2.1, 0.8, 1.0, 0.5, 0.5, 0.3, 0.2, 0.1)
results = fitModel(dataX = list(trials), dataResponse = list(responses))
fitParams = results$par[1,]
fitEst = computeModel(trials, mFast=fitParams[1], mSlow=fitParams[2], 
					  n=fitParams[3], h=fitParams[4], g=fitParams[7], 
					  verbose = TRUE)
plot(fitEst, actual = responses / max(responses))


###################################################
### code chunk number 4: response
###################################################
fitResponses = calculateResponse(fitParams[5], fitParams[6], fitEst)
plot(1:30, responses, type = "b", bty = "l",
	 xlab = "Trial", ylab = "Response rate")
points(1:30, fitResponses, col = 4, pch = 19)


