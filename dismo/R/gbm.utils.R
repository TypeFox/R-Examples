

.roc <-function (obsdat, preddat) {
# code adapted from Ferrier, Pearce and Watson's code, by J.Elith
#
# see:
# Hanley, J.A. & McNeil, B.J. (1982) The meaning and use of the area
# under a Receiver Operating Characteristic (ROC) curve.
# Radiology, 143, 29-36
#
# Pearce, J. & Ferrier, S. (2000) Evaluating the predictive performance
# of habitat models developed using logistic regression.
# Ecological Modelling, 133, 225-245.
# this is the non-parametric calculation for area under the ROC curve, 
# using the fact that a MannWhitney U statistic is closely related to
# the area
#

# in dismo, this is used in the gbm routines, but not elsewhere (see evaluate).

    if (length(obsdat) != length(preddat)) { 
        stop("obs and preds must be equal lengths")
	}
    n.x <- length(obsdat[obsdat == 0])
    n.y <- length(obsdat[obsdat == 1])
    xy <- c(preddat[obsdat == 0], preddat[obsdat == 1])
    rnk <- rank(xy)
    wilc <- ((n.x * n.y) + ((n.x * (n.x + 1))/2) - sum(rnk[1:n.x]))/(n.x * n.y)
    return(round(wilc, 4))
}




.calibration <- function(obs, preds, family = "binomial")  {
#
# j elith/j leathwick 17th March 2005
# calculates calibration statistics for either binomial or count data
# but the family argument must be specified for the latter 
# a conditional test for the latter will catch most failures to specify
# the family
#

	if (family == "bernoulli") {
		family <- "binomial"
	}
	
	pred.range <- max(preds) - min(preds)
	
	if(pred.range > 1.2 & family == "binomial") {
		print(paste("range of response variable is ", round(pred.range, 2)), sep = "", quote = F)
		print("check family specification", quote = F)
		return()
	}
	
	if(family == "binomial") {
		pred <- preds + 1e-005
		pred[pred >= 1] <- 0.99999
		mod <- glm(obs ~ log((pred)/(1 - (pred))), family = binomial)
		lp <- log((pred)/(1 - (pred)))
		a0b1 <- glm(obs ~ offset(lp) - 1, family = binomial)
		miller1 <- 1 - pchisq(a0b1$deviance - mod$deviance, 2)
		ab1 <- glm(obs ~ offset(lp), family = binomial)
		miller2 <- 1 - pchisq(a0b1$deviance - ab1$deviance, 1)
		miller3 <- 1 - pchisq(ab1$deviance - mod$deviance, 1)
	} else if(family == "poisson") {
		mod <- glm(obs ~ log(preds), family = poisson)
		lp <- log(preds)
		a0b1 <- glm(obs ~ offset(lp) - 1, family = poisson)
		miller1 <- 1 - pchisq(a0b1$deviance - mod$deviance, 2)
		ab1 <- glm(obs ~ offset(lp), family = poisson)
		miller2 <- 1 - pchisq(a0b1$deviance - ab1$deviance, 1)
		miller3 <- 1 - pchisq(ab1$deviance - mod$deviance, 1)
	}

	calibration.result <- c(mod$coef, miller1, miller2, miller3)
	names(calibration.result) <- c("intercept", "slope", "testa0b1", "testa0|b1", "testb1|a")
	return(calibration.result)
}

