library(discSurv)

#####################
# Unit Tests
#####################

######################
# gumbel

# Checks if appropriate names are specified in function
Check <- any(names(gumbel())=="linkfun") & 
   any(names(gumbel())=="linkinv") & 
   any(names(gumbel())=="mu.eta") & 
         any(names(gumbel())=="valideta") & 
               any(names(gumbel())=="name")
stopifnot(Check)
# Check if linkinv is an inverse function of linkfun by examples
ExpectGiven0 <- gumbel()$linkinv (0)
stopifnot(all.equal(0, gumbel()$linkfun (ExpectGiven0)))
ExpectGiven1 <- gumbel()$linkinv (1)
stopifnot(all.equal(1, gumbel()$linkfun (ExpectGiven1)))
ExpectGivenMinus1 <- gumbel()$linkinv (-1)
stopifnot(all.equal(-1, gumbel()$linkfun (ExpectGivenMinus1)))
# Check if numerical derivative and exact derivative are close together by examples
library(numDeriv)
stopifnot(all.equal(grad(func=gumbel()$linkinv, x=0), gumbel()$mu.eta (0)))
stopifnot(all.equal(grad(func=gumbel()$linkinv, x=1), gumbel()$mu.eta (1)))
stopifnot(all.equal(grad(func=gumbel()$linkinv, x=-1), gumbel()$mu.eta (-1)))

###########################
# estSurv

TrueProbs <- c(0.2, 0.3, 0.1, 0.3, 0.1)
TrueSurv <- sapply(1:length(TrueProbs), function (x) 1-sum(TrueProbs [1:x]))
TrueHaz <- c(TrueProbs [1], sapply(2:length(TrueProbs), function (x) TrueProbs [x] / TrueSurv [x-1]))
EstSurvFunc <- as.numeric(estSurv(TrueHaz [1:(length(TrueHaz)-1)]))
stopifnot(all.equal(EstSurvFunc, TrueSurv))

#############################
# estMargProb

TrueProbs <- (10:1)/sum(10:1)
TrueSurv <- sapply(1:length(TrueProbs), function (x) 1-sum(TrueProbs [1:x]))
TrueHaz <- c(TrueProbs [1], sapply(2:length(TrueProbs), function (x) TrueProbs [x] / TrueSurv [x-1]))
EstMargProb <- as.numeric(estMargProb(TrueHaz [1:(length(TrueHaz)-1)]))
stopifnot(all.equal(EstMargProb, TrueProbs))

# Case with only a single observed time interval
TrueProbs <- c(1/5, 4/5)
TrueSurv <- sapply(1:length(TrueProbs), function (x) 1-sum(TrueProbs [1:x]))
TrueHaz <- c(TrueProbs [1], sapply(2:length(TrueProbs), function (x) TrueProbs [x] / TrueSurv [x-1]))
EstMargProb <- as.numeric(estMargProb(TrueHaz [1:(length(TrueHaz)-1)]))
stopifnot(all.equal(EstMargProb, TrueProbs))