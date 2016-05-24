estate <-
function(Y,Z,X=NULL,Ypre=NULL,prob=NULL,HT=FALSE) {

# Raj difference if a pretest is available
	Yd <- Y
	if(!is.null(Ypre)) Yd <- Y - Ypre
	
	if(is.null(prob)) {
		warning("Probabilities not specified. Assuming equal probabilities.")
		prob <- mean(Z)
	}
	
	# does not work?

	weights <- Z/prob + (1-Z)/(1-prob)

# IPW
	if (HT == FALSE & is.null(X)) {
		ate <- weighted.mean(Yd[Z==1],weights[Z==1]) - weighted.mean(Yd[Z==0],weights[Z==0])

	}

# IPW + Covariates
	if (HT == FALSE & !is.null(X)) {
		ate <- lm(Yd~Z+X,weights = weights)$coefficients[2]
	}
	
# HT
	if (HT == TRUE & is.null(X)) ate <- (sum(Yd[Z==1]*weights[Z==1]) - sum(Yd[Z==0]*weights[Z==0]))/length(Yd)

# Raj-Within
	if (HT == TRUE & !is.null(X)) {
		Ypre2 <- lm(Yd~X)$fitted
		Yd2 <- Yd - Ypre2
		ate <- (sum(Yd2[Z==1]*weights[Z==1]) - sum(Yd2[Z==0]*weights[Z==0]))/length(Yd)
		}
	
	return(ate)
	}
