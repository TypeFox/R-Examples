# these functions originally by Jarrett Byrnes

# methods to generate 
# various information criteria from
# sem or adjchisq objects
# as well as generate and AIC table

# last modified 2012-10-02 by J. Fox

logLik.objectiveML <- function(object, ...){
	-0.5*deviance(object)
	}

# generics

AICc <- function(object, ...) UseMethod("AICc")

CAIC <- function(object, ...) UseMethod ("CAIC")


# methods for sem objects

AIC.objectiveML <- function(object, ..., k) {
	deviance(object) + 2*object$t
}

# small sample second order corrected aic
AICc.objectiveML <- function(object, ...) {
	deviance(object) + 2*object$t*(object$t + 1)/(object$N - object$t - 1)
}

# Consistent Akaike Information Criterion
CAIC.objectiveML <- function(object, ...) {
	props <- semProps(object)
	props$chisq - props$df*(1 + log(object$N))
}

BIC.objectiveML <- function(object, ...) {
    n <- object$n
    n.fix <- object$n.fix
    N <- object$N
    t <- object$t
    df <- n*(n + 1)/2 - t - n.fix*(n.fix + 1)/2
#	deviance(object) + object$t*log(object$N)
    deviance(object) - df*log(N)
}

# the following are not exported and are just place-keepers for the summary method

BIC.objectiveGLS <- function(object, ...) NULL

AIC.objectiveGLS <- function(object, ...) NULL

AICc.objectiveGLS <- function(object, ...) NULL

CAIC.objectiveGLS <- function(object, ...) NULL

# weights

aicW <- function(a.list, func=AICc){
	aiclist <- sapply(a.list, function(x) eval(func(x)), simplify=TRUE)			
	delta.i <- aiclist - min(aiclist)
	aicw <- exp(-0.5*delta.i)/sum(exp(-0.5*delta.i))
	return.matrix <- matrix(c(aiclist,delta.i, aicw), ncol=3)
	colnames(return.matrix) <- c("IC", "delta.i", "weight")
	rownames(return.matrix) <- 1:length(return.matrix[,1])
	return(return.matrix)
}
