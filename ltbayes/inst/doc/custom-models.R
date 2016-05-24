## ----custom3pl-----------------------------------------------------------
custom3pl <- function(zeta, y, alph, beta, gamm, prior = dnorm, ...) {
	if (is.vector(y)) y <- matrix(y, 1, length(y))
	m <- ncol(y)
	n <- nrow(y)
	prob <- matrix(NA, m, 2)
	prob[,2] <- gamm + (1 - gamm) * plogis(alph * (zeta - beta))
	prob[,1] <- 1 - prob[,2]
	yprb <- matrix(NA, n, m)
	for (i in 1:n) {
		yprb[i,] <- prob[row(prob) == 1:m & col(prob) == y[i,] + 1]
	}
	return(list(post = log(sum(apply(yprb, 1, prod))) 
		+ log(prior(zeta, ...)), prob = prob))
}

## ----, eval = FALSE------------------------------------------------------
#  custom3pl <- function(zeta, y, alph, beta, prior = dnorm, ...) {
#  	prob <- gamm + (1 - gamm) * plogis(alph * (zeta - beta))
#  	yprb <- prob^y * (1 - prb)^(1 - y)
#  	return(list(post = log(prob(yprb))
#  		+ log(prior(zeta, ...)), prob = prob))
#  }	

## ----, echo = FALSE, message = FALSE-------------------------------------
library(ltbayes) 

## ----fmodel3pldemo-------------------------------------------------------
samp <- 5000 # number of sampled realizations from posterior distribution
burn <- 1000 # number of discarded burn-in samples
alph <- c(1.27,1.34,1.14,1,0.67)   # discrimination parameters
beta <- c(1.19,0.59,0.15,-0.59,-2) # difficulty parameters
gamm <- c(0.1,0.15,0.15,0.2,0.1)   # guessing parameters
set.seed(123)
zeta.fmodel3pl <- postsamp(fmodel3pl, c(0,0,1,1,1), 
	apar = alph, bpar = beta, cpar = gamm,
	control = list(nbatch = samp + burn, scale = 3))
zeta.fmodel3pl <- data.frame(sample = 1:samp, 
	zeta = zeta.fmodel3pl$batch[(burn + 1):(samp + burn)])
head(zeta.fmodel3pl)
set.seed(123)
zeta.custom3pl <- postsamp(fmodel3pl, c(0,0,1,1,1), 
	apar = alph, bpar = beta, cpar = gamm,
	control = list(nbatch = samp + burn, scale = 3))
zeta.custom3pl <- data.frame(sample = 1:samp, 
	zeta = zeta.custom3pl$batch[(burn + 1):(samp + burn)])
head(zeta.custom3pl)

## ----, gpc---------------------------------------------------------------
fmodelgpc <- function(zeta, y, apar, bpar, prior = dnorm, ...) {
	if (is.vector(y)) y <- matrix(y, 1, length(y))
	m <- ncol(y)
	n <- nrow(y)
	r <- ncol(beta) + 1
	prob <- exp(outer(apar*zeta, 0:(r-1)) - 
		t(apply(sweep(cbind(0, bpar), 1, apar, "*"), 1, cumsum)))
	prob <- sweep(prob, 1, apply(prob, 1, sum), "/")	
	yprb <- matrix(NA, n, m)
	for (i in 1:n) {
		yprb[i,] <- prob[row(prob) == 1:m & col(prob) == y[i,] + 1]
	}
	return(list(post = log(sum(apply(yprb, 1, prod)))
		+ log(prior(zeta, ...)), prob = prob))
}

## ----gpcdemo-------------------------------------------------------------
alph <- rep(1, 5)        # "discrimination" parameters
beta <- matrix(0, 5, 2)  # "difficulty" parameters
set.seed(123)
zeta.fmodelpcm <- postsamp(fmodelpcm, c(0,1,2,1,0), bpar = beta,
	control = list(nbatch = samp + burn))
zeta.fmodelpcm <- data.frame(sample = 1:samp, 
	zeta = zeta.fmodelpcm$batch[(burn + 1):(samp + burn)])
head(zeta.fmodelpcm)
set.seed(123)
zeta.fmodelgpc <- postsamp(fmodelgpc, c(0,1,2,1,0), apar = alph, bpar = beta,
	control = list(nbatch = samp + burn))
zeta.fmodelgpc <- data.frame(sample = 1:samp, 
	zeta = zeta.fmodelgpc$batch[(burn + 1):(samp + burn)])
head(zeta.fmodelgpc)

