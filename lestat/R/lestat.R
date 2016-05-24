
######################
### Generic functions: 
######################

expectation <- function(object) {
UseMethod("expectation")
}

expectation.default <- function(object) {
print("The expectation function is not implemented for this object.")
}

variance <- function(object) {
UseMethod("variance")
}

variance.default <- function(object) {
print("The variance function is not implemented for this object.")
}

precision <- function(object) {
UseMethod("precision")
}

precision.default <- function(object) {
cat("The precision function is not implemented for this object.\n")
}

marginal <- function(object, v) {
UseMethod("marginal")
}

marginal.default <- function(object, v) {
print("The marginal function is not implemented for this object.")
}

conditional <- function(object, v, val) {
UseMethod("conditional")
}

conditional.default <- function(object, v, val) {
print("The conditional function is not implemented for this object.")
}

probability <- function(object, val) {
UseMethod("probability")
}

probability.default <- function(object, val) {
print("The probability function is not implemented for this object.")
}

probabilitydensity <- function(object, val, log=FALSE, normalize=TRUE) {
UseMethod("probabilitydensity")
}

probabilitydensity.default <- function(object, val, log=FALSE, normalize=TRUE) {
print("The probabilitydensity function is not implemented for this object.")
}

cdf <- function(object, val) {
UseMethod("cdf")
}

cdf.default <- function(object, val) {
cat("The cdf function is not implemented for this object.\n")
}

invcdf <- function(object, val) {
UseMethod("invcdf")
}

invcdf.default <- function(object, val) {
cat("The invcdf function is not implemented for this object.\n")
}

# ALREADY GENERIC: 
# print(x,...)
# plot(x, y, ...)
# summary(object, ...)
# simulate(object, nsim = 1, seed = NULL, ...)


difference <- function(object1, object2) {
# This function assumes that object1 and object2 are independent distributions, 
# and if they are normal distributions or t-distributions, generates the (approximate)
# distribution representing the difference. 
# This should mainly be used to generate the APPROXIMATE t-distribution resulting from 
# taking the difference of two 1D t-distributions. 
UseMethod("difference")
}

difference.default <- function(object1, object2) {
print("The difference function is not implemented for this combination of objects.")
}

compose <- function(object, type, ...) {
   UseMethod("compose")
}

compose.default <- function(object, type, ...) {
   cat("The compose function is not implemented for this object.\n")
}

linearpredict <- function(object, ...) {
UseMethod("linearpredict")
}

linearpredict.default <- function(object, ...) {
   cat("The linearpredict function is not implemented for this object.\n")
}

prediction <- function(object) {
UseMethod("prediction")
}

prediction.default <- function(object) {
   cat("The prediction function is not implemented for this object.\n")
}

contrast <- function(object, v) {
UseMethod("contrast")
}

contrast.default <- function(object, v) {
   cat("The contrast function is not implemented for this object.\n")
}

credibilityinterval <- function(object, prob=0.95) {
UseMethod("credibilityinterval")
}

credibilityinterval.default <- function(object, prob) {
cat("The credibilityinterval function is not implemented for this object.\n")
}

p.value <- function(object, point) {
# For any distribution, this should return the probability outside the 
# credibility region with boundary containing the origin, or if value or vector
# is given, that value or vector.  
UseMethod("p.value")
}

p.value.default <- function(object, point) {
print("The p.value function is not implemented for this object.")
}

##################################
## CLASS discretedistribution 
##################################
# Data: 
# vals
# probs

discretedistribution <- function(vals, probs = rep(1, length(vals))) {
   if (length(unique(vals))<length(vals)) 
      cat("ERROR: The discrete distribution cannot be made with repeated elements.\n")
   else if (length(probs)!=length(vals)) 
      cat("ERROR: Wrong length for the probability vector.\n")
   else if (min(probs)<=0)
      cat("ERROR: All probabilities must be positive.\n")
   else {
      result <- list(vals = vals, probs = probs/sum(probs))
      class(result) <- c("discretedistribution", "probabilitydistribution")
      result
   }
}

expectation.discretedistribution <- function(object) {
   if (is.numeric(object$vals)) 
      sum(object$vals*object$probs)
   else 
      cat("ERROR: Cannot compute the expectation when outcomes are non-numeric.\n")
}

variance.discretedistribution <- function(object) {
   if (is.numeric(object$vals)) 
      sum(object$vals^2*object$probs) - sum(object$vals*object$probs)^2
   else 
      cat("ERROR: Cannot compute the variance when outcomes are non-numeric.\n")
}

precision.discretedistribution <- function(object) {
   if (is.numeric(object$vals)) 
      1/(sum(object$vals^2*object$probs) - sum(object$vals*object$probs)^2)
   else 
      cat("ERROR: Cannot compute the precision when outcomes are non-numeric.\n")
}

probability.discretedistribution <- function(object, val) {
   index <- val==object$vals
   if (sum(index)==1) 
      object$probs[index]
   else
      cat("ERROR: Given value not a possible value for the distribution.\n")
}

cdf.discretedistribution <- function(object, val) {
   index <- val==object$vals
   if (sum(index)==1) 
      sum(object$probs[1:((1:length(object$vals))[index])])
   else
      cat("ERROR: Given value not a possible value for the distribution.\n")
}

invcdf.discretedistribution <- function(object, val) {
   if (val < 0) 
      cat("ERROR: Second argument cannot be negative.\n")
   else if (val>1)
      cat("ERROR: Second argument cannto be greater than 1.\n")
   else 
       object$vals[sum(cumsum(object$probs)<=val)+1]
}

print.discretedistribution <- function(x, ...) {
   cat("A discrete probability distribution.\n")
   names(x$probs) <- x$vals
   print(x$probs)
}

plot.discretedistribution <- function(x, ...) {
   n <- length(x$vals)
   plot(x$probs, type="h", ylab="Probability")
   text(labels=x$vals, x=1:n + 0.01, y=x$probs)
}

summary.discretedistribution <- function(object, ...) {
   print(object)
}

simulate.discretedistribution <- function(object, nsim = 1, ...) {
   sample(object$vals, nsim, replace = TRUE, prob = object$probs)
}

compose.discretedistribution <- function(object, type, ...) {
   if (type!="discretedistribution") 
      cat("ERROR: Not implemented for this type.\n")
   else {
      args <- list(...) 
      mdiscretedistribution(diag(object$probs)%*%args[[2]], list(X1 = object$vals, X2 = args[[1]]))
   }
}

##################################
## CLASS mdiscretedistribution 
##################################

mdiscretedistribution <- function(probs, nms=NULL) {
   if (is.null(dim(probs)))
      cat("ERROR: The probabilities must be in a matrix or an array.\n")
   else {
      if (is.null(nms)) {
      	 nms <- list()
         for (i in 1:length(dim(probs))) nms[[i]] <- 1:(dim(probs)[i])
      }
      dimnames(probs) <- nms
      result <- list(vals = nms, probs = probs/sum(probs))
      class(result) <- c("mdiscretedistribution", "probabilitydistribution")
      result
   }
}

expectation.mdiscretedistribution <- function(object) {
   if (!is.numeric(unlist(object$vals))) 
      cat("ERROR: Cannot compute the expectation of this distribution.\n")
   else {
      k <- length(dim(object$probs))
      result <- rep(NA, k)
      for (i in 1:k) result[i] <- expectation(marginal(object, i))
      result
   } 
}

variance.mdiscretedistribution <- function(object) {
   if (is.numeric(unlist(object$vals))) {
      k <- length(dim(object$probs))
      result <- matrix(NA, k, k)
      for (i in 1:k) 
         for (j in 1:k) {
	    if (i==j) 
	       result[i,j] <- variance(marginal(object, i))
            else {
	       XX <- outer(object$vals[[i]], object$vals[[j]])
               YY <- apply(object$probs, c(i,j), sum)
	       result[i,j] <- sum(XX*YY) - expectation(marginal(object, i))*expectation(marginal(object, j))
            }
      }
      result
   } else     
      cat("ERROR: Cannot compute the variance when outcomes are non-numeric.\n")
}

precision.mdiscretedistribution <- function(object) {
   if (is.numeric(unlist(object$vals)))
      solve(variance(object))
   else 
      cat("ERROR: Cannot compute the precision when outcomes are non-numeric.\n")
}

probability.mdiscretedistribution <- function(object, val) {
   k <- length(dim(object$probs))
   index <- 0
   for (i in k:1) 
      index  <- index*dim(object$probs)[i] + (1:dim(object$probs)[i])[object$vals[[i]]==val[i]] - 1
   object$probs[index+1]
}

print.mdiscretedistribution <- function(x, ...) {
   cat("A multivariate discrete probability distribution.\n")
   print(x$probs)
}

plot.mdiscretedistribution <- function(x, onlybivariate=FALSE, ...) {
   if (onlybivariate & length(dim(x$probs))==2) {
      xx <- x$vals[[1]]
      if (!is.numeric(xx)) xx <- 1:(length(xx))
      yy <- x$vals[[2]]
      if (!is.numeric(yy)) yy <- 1:(length(yy))
      image(xx, yy, x$probs, xlab="", ylab="")
   } else {
      k <- length(dim(x$probs))
      old.par <- par(no.readonly = TRUE); on.exit(par(old.par))
      par(mfcol=c(k,k))
      for (i in 1:k) for (j in 1:k) 
         if (i==j) plot(marginal(x, i))
         else plot(marginal(x, c(i,j)), TRUE)
   }
}

summary.mdiscretedistribution <- function(object, ...) {
   print(object)
}

simulate.mdiscretedistribution <- function(object, nsim = 1, ...) {
   n <- length(object$probs)
   k <- length(dim(object$probs))
   ind <- sample(0:(n-1), nsim, replace = TRUE, prob = object$probs)
   result <- matrix(NA, nsim, k)
   for (i in 1:k) {
      result[,i] <- object$vals[[i]][ind %% dim(object$probs)[i] + 1]
      ind <- ind %/% dim(object$probs)[i] 
   }
   result
}

marginal.mdiscretedistribution <- function(object, v)
{
   object$probs <- apply(object$probs, v, sum)
   if (length(v)==1) 
      discretedistribution(object$vals[[v]], object$probs)
   else
   {  
      resultvals <- list()
      for (i in 1:length(v)) resultvals[[i]] <- object$vals[[v[i]]]
      object$vals <- resultvals
      object
   }
}

conditional.mdiscretedistribution <- function(object, v, val) 
{
   n <- length(object$probs)
   k <- length(dim(object$probs))
   index <- rep(TRUE, n)
   counter <- 1
   for (i in v) {
      ind <- (1:dim(object$probs)[i])[object$vals[[i]]==val[counter]] - 1
      ind2 <- 0:(n-1)
      if (i>1) ind2 <- ind2%/%prod(dim(object$probs)[1:(i-1)])
      ind2 <- ind2%%dim(object$probs)[i]
      index <- index & (ind2==ind)
      counter <- counter + 1
   }
   result <- object$probs[index]
   dim(result) <- dim(object$probs)[-v]
   resultvals <- object$vals
   v <- sort(v, decreasing=TRUE)
   for (i in v) 
      resultvals[[i]] <- NULL
   if (k-length(v)==1) 
      discretedistribution(resultvals[[1]], result)
   else 
      mdiscretedistribution(result, resultvals)
}

##################################
## CLASS binomialdistribution 
##################################
# Data: 
# count, probability
# OPTIONAL: name: of the single dimension

binomialdistribution <- function(ntrials, probability) {
   if (ntrials <1) 
      cat("ERROR: The first parameter must be a positive integer.\n")
   else if (probability < 0) 
      cat("ERROR: The probability cannot be smaller than zero.\n")
   else if (probability > 1) 
      cat("ERROR: The probability cannot be larger than 1.\n")
   else {
      result <- list(ntrials=ntrials, probability=probability)
      class(result) <- c("binomialdistribution", "probabilitydistribution")
      result
   }
}

expectation.binomialdistribution <- function(object) {
   object$ntrials*object$probability
}

variance.binomialdistribution <- function(object) {
   object$ntrials*object$probability*(1-object$probability)
}

precision.binomialdistribution <- function(object) {
   1/(object$ntrials*object$probability*(1-object$probability))
}

probability.binomialdistribution <- function(object, val) {
   dbinom(val, object$ntrials, object$probability)
}

cdf.binomialdistribution <- function(object, val) {
   pbinom(val, object$ntrials, object$probability)
}

invcdf.binomialdistribution <- function(object, val) {
   qbinom(val, object$ntrials, object$probability)
}

print.binomialdistribution <- function(x, ...) {
   cat("A Binomial probability distribution.\n")
   cat("Number of trials ", x$ntrials, ", probability ", x$probability, ".\n", sep="")
}

plot.binomialdistribution <- function(x, ...) {
   plot(0:x$ntrials, dbinom(0:x$ntrials, x$ntrials, x$probability), xlab="", ylab="")
}

summary.binomialdistribution <- function(object, ...) {
   cat("A Binomial probability distribution.\n")
   cat("Number of trials ", object$ntrials, ", probability ", object$probability, ".\n", sep="")
}

simulate.binomialdistribution <- function(object, nsim = 1, ...) {
   rbinom(nsim, object$ntrials, object$probability)
}

##################################
## CLASS betadistribution
##################################

betadistribution <- function(alpha, beta) {
      result <- list(alpha=alpha, beta=beta)
      class(result) <- c("betadistribution", "probabilitydistribution")
      result
}

expectation.betadistribution <- function(object) {
   object$alpha/(object$alpha + object$beta)
}

variance.betadistribution <- function(object) {
   object$alpha*object$beta/((object$alpha + object$beta)^2*(object$alpha + object$beta + 1))
}

precision.betadistribution <- function(object) {
   1/variance(object)
}

probabilitydensity.betadistribution <- function(object, val, log=FALSE, normalize=TRUE) {
   dbeta(val, object$alpha, object$beta)
}

cdf.betadistribution <- function(object, val) {
   pbeta(val, object$alpha, object$beta)
}

invcdf.betadistribution <- function(object, val) {
   qbinom(val, object$alpha, object$beta)
}

print.betadistribution <- function(x, ...) {
   cat("A Beta probability distribution.\n")
   cat("Alpha ", x$alpha, ", beta ", x$beta, ".\n", sep="")
}

plot.betadistribution <- function(x, ...) {
   xx <- (1:999)/1000
   plot(xx, dbeta(xx, x$alpha, x$beta), type="l", xlab="", ylab="")
}

summary.betadistribution <- function(object, ...) {
   cat("A Beta probability distribution.\n")
   cat("Alpha ", object$alpha, ", beta ", object$beta, ".\n", sep="")
}

simulate.betadistribution <- function(object, nsim = 1, ...) {
   rbeta(nsim, object$alpha, object$beta)
}

credibilityinterval.betadistribution <- function (object, prob = 0.95) 
{
    c(qbeta((1 - prob)/2, object$alpha, object$beta), 
        qbeta(1 - (1 - prob)/2, object$alpha, object$beta))
}

##################################
## CLASS multiomialdistribution 
##################################
# Data: 
# count, probability
# OPTIONAL: name: of the dimensions

multinomialdistribution <- function(ntrials, probabilities) {
   if (ntrials <1) 
      cat("ERROR: The first parameter must be a positive integer.\n")
   else if (any(probabilities < 0)) 
      cat("ERROR: The probabilities cannot be smaller than zero.\n")
   else if (sum(probabilities)==0)
      cat("ERROR: The probabilities cannot all be zero.\n")
   else {
      result <- list(ntrials=ntrials, probabilities = probabilities/sum(probabilities))
      class(result) <- c("multinomialdistribution", "probabilitydistribution")
      result
   }
}

expectation.multinomialdistribution <- function(object) {
   object$ntrials*object$probabilities
}

variance.multinomialdistribution <- function(object) {
   object$ntrials*(diag(object$probabilities) - outer(object$probabilities, object$probabilities))
}

probability.multinomialdistribution <- function(object, val) {
   if (sum(val)==object$ntrials) 
       dmultinom(val, object$ntrials, object$probabilities)
   else 0
}

print.multinomialdistribution <- function(x, ...) {
   cat("A Multinomial probability distribution.\n")
   cat("Number of trials ", x$ntrials, ", probabilities ", x$probabilities, ".\n", sep="")
}

#plot.multinomialdistribution <- function(x, ...) {
#   plot(0:x$ntrials, dbinom(0:x$ntrials, x$ntrials, x$probability), xlab="", ylab="")
#}

summary.multinomialdistribution <- function(object, ...) {
   cat("A Multinomial probability distribution.\n")
   cat("Number of trials ", object$ntrials, ", probabilities ", object$probabilities, ".\n", sep="")
}

simulate.multinomialdistribution <- function(object, nsim = 1, ...) {
   rmultinom(nsim, object$ntrials, object$probability)
}

##################################
## CLASS betadistribution
##################################

betadistribution <- function(alpha, beta) {
      result <- list(alpha=alpha, beta=beta)
      class(result) <- c("betadistribution", "probabilitydistribution")
      result
}

expectation.betadistribution <- function(object) {
   object$alpha/(object$alpha + object$beta)
}

variance.betadistribution <- function(object) {
   object$alpha*object$beta/((object$alpha + object$beta)^2*(object$alpha + object$beta + 1))
}

precision.betadistribution <- function(object) {
   1/variance(object)
}

probabilitydensity.betadistribution <- function(object, val, log=FALSE, normalize=TRUE) {
   dbeta(val, object$alpha, object$beta)
}

cdf.betadistribution <- function(object, val) {
   pbeta(val, object$alpha, object$beta)
}

invcdf.betadistribution <- function(object, val) {
   qbinom(val, object$alpha, object$beta)
}

print.betadistribution <- function(x, ...) {
   cat("A Beta probability distribution.\n")
   cat("Alpha ", x$alpha, ", beta ", x$beta, ".\n", sep="")
}

plot.betadistribution <- function(x, ...) {
   xx <- (1:999)/1000
   plot(xx, dbeta(xx, x$alpha, x$beta), type="l", xlab="", ylab="")
}

summary.betadistribution <- function(object, ...) {
   cat("A Beta probability distribution.\n")
   cat("Alpha ", object$alpha, ", beta ", object$beta, ".\n", sep="")
}

simulate.betadistribution <- function(object, nsim = 1, ...) {
   rbeta(nsim, object$alpha, object$beta)
}


##################################
## CLASS betabinomial
##################################
# this is the univariate marginal of the binomialbeta

betabinomial <- function(n, alpha, beta) {
      result <- list(n=n, alpha=alpha, beta=beta)
      class(result) <- c("betabinomial", "probabilitydistribution")
      result
}

expectation.betabinomial <- function(object) {
   object$n*object$alpha/(object$alpha + object$beta)
}

variance.betabinomial <- function(object) {
   object$n*object$alpha*object$beta*(object$alpha+object$beta+object$n)/
   ((object$alpha + object$beta)^2*(object$alpha + object$beta + 1))
}

precision.betabinomial <- function(object) {
   1/variance(object)
}

probability.betabinomial <- function(object, val) {
   gamma(object$n+1)*gamma(object$alpha+val)*gamma(object$n+object$beta-val)*gamma(object$alpha+object$beta)/
   (gamma(val+1)*gamma(object$n+1-val)*gamma(object$alpha+object$beta+object$n)*gamma(object$alpha)*gamma(object$beta))
}

cdf.betabinomial <- function(object, val) {
   sum(probability(object, 0:val))
}

#invcdf.betabinomial <- function(object, val) {
#   qbinom(val, object$alpha, object$beta)
#}

print.betabinomial <- function(x, ...) {
   cat("A Betabinomial probability distribution.\n")
   cat("n ", x$n, ", alpha ", x$alpha, ", beta ", x$beta, ".\n", sep="")
}

plot.betabinomial <- function(x, ...) {
   plot(0:x$n, probability(x, 0:x$n), xlab="", ylab="")
}

summary.betabinomial <- function(object, ...) {
   cat("A Betabinomial probability distribution.\n")
   cat("n ", object$n, ", alpha ", object$alpha, ", beta ", object$beta, ".\n", sep="")
}

simulate.betabinomial <- function(object, nsim = 1, ...) {
   simulate(binomialbeta(object$n, object$alpha, object$beta), nsim)[,2]
}

##################################
## CLASS binomialbeta
##################################
# this is the bivariate distribution combining the binomial with the beta

binomialbeta <- function(n, alpha, beta) {
      result <- list(n=n, alpha=alpha, beta=beta)
      class(result) <- c("binomialbeta", "probabilitydistribution")
      result
}

#expectation.binomialbeta <- function(object) {
#   object$n*object$alpha/(object$alpha + object$beta)
#}

#variance.binomialbeta <- function(object) {
#   object$n*object$alpha*object$beta*(object$alpha+object$beta+object$n)/
#   ((object$alpha + object$beta)^2*(object$alpha + object$beta + 1))
#}

#precision.binomialbeta <- function(object) {
#   1/variance(object)
#}

#probability.binomialbeta <- function(object, val) {
#   dbeta(val, object$alpha, object$beta)
#}

print.binomialbeta <- function(x, ...) {
   cat("A Binomialbeta probability distribution.\n")
   cat("n ", x$n, ", alpha ", x$alpha, ", beta ", x$beta, ".\n", sep="")
}

plot.binomialbeta <- function(x, onlybivariate=FALSE, switchaxes=FALSE, ...) {
# The onlybivariate parameter can make the output consist of only the bivariate plot
# The switchaxes parameter can make the axes switch, if the plot is only bivariate. 
   if (onlybivariate) {
      xx <- (1:999)/1000
      betad <- dbeta(rep(xx, x$n + 1), x$alpha, x$beta)
      yy <- 0:x$n
      binomd <- dbinom(rep(0:x$n, each=999), x$n, rep(xx, x$n + 1))
      zz <- matrix(betad*binomd, 999, x$n + 1 )
      if (switchaxes) 
      	 image(xx, yy, zz, xlab="", ylab="")
      else
	 image(yy ,xx, t(zz), xlab="", ylab="")
   } else {
      old.par <- par(no.readonly = TRUE); on.exit(par(old.par))
      par(mfcol=c(2,2))
      plot(marginal(x, 1))
      plot(x, TRUE, TRUE)
      plot(x, TRUE, FALSE)
      plot(marginal(x, 2))
   }
}

summary.binomialbeta <- function(object, ...) {
   cat("A Binomialbeta probability distribution.\n")
   cat("n ", object$n, ", alpha ", object$alpha, ", beta ", object$beta, ".\n", sep="")
}

simulate.binomialbeta <- function(object, nsim = 1, ...) {
   x <- rbeta(nsim, object$alpha, object$beta)
   y <- rbinom(nsim, object$n, x)
   cbind(x,y)
}

marginal.binomialbeta <- function(object, v) 
{
if (v==1) 
   betadistribution(object$alpha, object$beta)
else if (v==2)
   betabinomial(object$n, object$alpha, object$beta)
else
   cat("ERROR: Wrong specification.\n")
}

conditional.binomialbeta <- function(object, v, val)
{
if (v==1) 
   binomialdistribution(object$n, val)
else if (v==2)
   betadistribution(object$alpha + val, object$beta + object$n - val)
else 
   cat("ERROR: Wrong specification.\n")
}

##################################
## CLASS poissondistribution 
##################################
# Data: 
# rate
# OPTIONAL: name: of the single dimension

poissondistribution <- function(rate) {
   if (rate < 0) 
      cat("ERROR: The rate cannot be smaller than zero.\n")
   else {
      result <- list(rate=rate)
      class(result) <- c("poissondistribution", "probabilitydistribution")
      result
   }
}

expectation.poissondistribution <- function(object) {
   object$rate
}

variance.poissondistribution <- function(object) {
   object$rate
}

precision.poissondistribution <- function(object) {
   1/object$rate
}

probability.poissondistribution <- function(object, val) {
   dpois(val, object$rate)
}

cdf.poissondistribution <- function(object, val) {
   ppois(val, object$rate)
}

invcdf.poissondistribution <- function(object, val) {
   qpois(val, object$rate)
}

print.poissondistribution <- function(x, ...) {
   cat("A Poisson probability distribution.\n")
   cat("Rate ", x$rate, ".\n", sep="")
}

plot.poissondistribution <- function(x, ...) {
   plot(0:(x$rate*4), dpois(0:(x$rate*4), x$rate), xlab="", ylab="")
}

summary.poissondistribution <- function(object, ...) {
   cat("A Poisson probability distribution.\n")
   cat("Rate ", object$rate, ".\n", sep="")
}

simulate.poissondistribution <- function(object, nsim = 1, ...) {
   rpois(nsim, object$rate)
}


##################################
## CLASS uniformdistribution 
##################################
# Data: 
# a, b. 

uniformdistribution <- function(a=0, b=1) {
   if (a >= b) 
      cat("ERROR: The lower limit must be smaller than the upper limit.\n")
   else {
      result <- list(a=a, b=b)
      class(result) <- c("uniformdistribution", "probabilitydistribution")
      result
   }
}

expectation.uniformdistribution <- function(object) {
   (object$a+object$b)/2
}

variance.uniformdistribution <- function(object) {
   (object$b-object$a)^2/12
}

precision.uniformdistribution <- function(object) {
   12/(object$b-object$a)^2
}

probabilitydensity.uniformdistribution <- function(object, val, log=FALSE, normalize=TRUE) {
   if (val < object$a | val > object$b)
      ifelse(log, log(0), 0) 
   else if (normalize)
      ifelse(log, -log(object$b-object$a), 1/(object$b-object$a))
   else 
      ifelse(log, 0, 1)
}

cdf.uniformdistribution <- function(object, val) {
   if (val < object$a) 
     0 
   else if (val > object$b) 
     1 
   else 
     (val-object$a)/(object$b-object$a)
}

invcdf.uniformdistribution <- function(object, val) {
   if (val <0 | val > 1) 
      cat("ERROR: Wrong specification\n")
   else
      val*(object$b-object$a) + object$a
}

print.uniformdistribution <- function(x, ...) {
   cat("A Uniform probability distribution\n")
   cat("on the interval from ", x$a, " to ", x$b, ".\n", sep="")
}

plot.uniformdistribution <- function(x, ...) {
   if ((x$b - x$a)==Inf)
      cat("The distribution is improper.\n")
   else {
      xmin <- x$a - 0.1*(x$b-x$a)
      xmax <- x$b + 0.1*(x$b-x$a)
      plot(c(xmin, x$a,x$a,x$b,x$b, xmax), 
      c(0, 0, 1/(x$b-x$a), 1/(x$b-x$a), 0, 0), type="l", xlab="", ylab="")
}
}

summary.uniformdistribution <- function(object, ...) {
   print(object)
}

simulate.uniformdistribution <- function(object, nsim = 1, ...) {
   runif(nsim, object$a, object$b)
}

compose.uniformdistribution <- function(object, type, ...) {
# NOTE: This must be re-written as the number of types expands...
   if (type!="binomialdistribution") 
      cat("ERROR: Not implemented for this type.\n")
   else if (object$a != 0 | object$b != 1)
      cat("ERROR: The uniform distribution must be on the interval from 0 to 1.\n")
   else {
      args <- list(...)
      binomialbeta(args[[1]], 1, 1)
   }
}

##################################
## CLASS muniformdistribution 
##################################

muniformdistribution <- function(startvec, stopvec) {
   if (any(startvec>=stopvec)) 
      cat("ERROR: The lower limit must be smaller than the upper limit.\n")
   else {
      result <- list(a=startvec, b=stopvec)
      class(result) <- c("muniformdistribution", "probabilitydistribution")
      result
   }
}

expectation.muniformdistribution <- function(object) {
   (object$a+object$b)/2
}

variance.muniformdistribution <- function(object) {
   diag((object$b-object$a)^2/12)
}

precision.muniformdistribution <- function(object) {
   diag(12/(object$b-object$a)^2)
}

probabilitydensity.muniformdistribution <- function(object, val, log=FALSE, normalize=TRUE) {
   if (any(val < object$a) | any(val > object$b))
      ifelse(log, log(0), 0) 
   else if (normalize)
      ifelse(log, -sum(log(object$b-object$a)), 1/prod(object$b-object$a))
   else 
      ifelse(log, 0, 1)
}

print.muniformdistribution <- function(x, ...) {
   cat("A Multivariate uniform probability distribution on intervals\n")
   outp <- cbind(x$a, x$b)
   colnames(outp) <- c("start", "stop")
   rownames(outp) <- NULL
   print(outp)
}

plot.muniformdistribution <- function(x, onlybivariate=FALSE, ...) {
   if (any((x$b - x$a)==Inf))
      cat("The distribution is improper.\n")
   else if (onlybivariate & length(x$a)==2) {
      xmin <- x$a - 0.1*(x$b-x$a)
      xmax <- x$b + 0.1*(x$b-x$a)
      plot(rbind(xmin, xmax), xlab="", ylab="", type="n")
      lines(matrix(c(x$a[1], x$a[1], x$b[1], x$b[1], x$a[1], x$a[2], x$b[2], x$b[2], x$a[2], x$a[2]), 5, 2))
   } else {
      k <- length(x$a)
      old.par <- par(no.readonly = TRUE); on.exit(par(old.par))
      par(mfcol=c(k,k))
      for (i in 1:k) for (j in 1:k) 
         if (i==j) plot(marginal(x, i))
         else plot(marginal(x, c(i,j)), TRUE)
   }
}

summary.muniformdistribution <- function(object, ...) {
   print(object)
}

simulate.muniformdistribution <- function(object, nsim = 1, ...) {
   matrix(runif(nsim*length(object$a), object$a, object$b), nsim, length(object$a), byrow=TRUE)
}

marginal.muniformdistribution <- function(object, v) { 
   # v gives the indices for the dimensions which should be kept. 
   if (length(v)==1)
      uniformdistribution(object$a[v], object$b[v])
   else
      muniformdistribution(object$a[v], object$b[v])
}

conditional.muniformdistribution <- function(object, v, val) {
# v gives the indices for the dimensions whose values should be fixed, and the 
# fixed values are given by val. 
   a <- object$a[-v]
   b <- object$b[-v]
   if (length(a)==1)
      uniformdistribution(a,b)
   else
      muniformdistribution(a,b)
}

##################################
## CLASS normal 
##################################
# Data: 
# expectation
# lambda: the logged scale, i.e., log of the standard deviation. 
# OPTIONAL alternative parameter: The precision P
# OPTIONAL: name: of the single dimension

normal <- function(expectation=0, lambda, P=1) {
   if (!missing(lambda) & !missing(P)) 
      cat("ERROR: You should not specify both logged scale and precision.\n")
   else if (!missing(P) && P<0) 
      cat("ERROR: The precision cannot be negative.\n")
   else {
      if (!missing(lambda)) P <- exp(-2*lambda)
      result <- list(expectation=expectation, P=P)
      class(result) <- c("normal", "probabilitydistribution")
      result
   }
}

expectation.normal <- function(object) {
   object$expectation
}

variance.normal <- function(object) {
   1/object$P
}

precision.normal <- function(object) {
   object$P
}

probabilitydensity.normal <- function(object, val, log=FALSE, normalize=TRUE) {
   if (normalize)
      dnorm(val, object$expectation, 1/sqrt(object$P), log)
   else
      dnorm(val, object$expectation, 1/sqrt(object$P), log)
}

cdf.normal <- function(object, val) {
   pnorm(val, object$expectation, 1/sqrt(object$P))
}

invcdf.normal <- function(object, val) {
   if (val<=0 | val>=1) 
      cat("ERROR: Illegal input for the invcdf function.\n")
   else 
      qnorm(val, object$expectation, 1/sqrt(object$P))
}

print.normal <- function(x, ...) {
if (x$P > 0)
{
cat("A univariate Normal probability distribution.\n")
cat("Expectation ", x$expectation, ", logged scale ", -0.5*log(x$P), ", SD ", sqrt(1/x$P), ".\n", sep="")
}
else
cat("A flat univariate Normal probability distribution.\n")
}

plot.normal <- function(x, ...) {
if (x$P > 0) {
   stdev <- sqrt(1/x$P)
   start <- qnorm(0.001, x$expectation, stdev)
   stop  <- qnorm(0.999, x$expectation, stdev)
   xx    <- seq(start, stop, length=1001)
   yy    <- dnorm(xx, mean = x$expectation, sd = stdev)
   plot(xx, yy, xlab = "", ylab="", type="l")
} else {
   plot(c(-3, 3), c(0, 2), xlab="", ylab="", type="n")
   abline(h=1)
}
}

summary.normal <- function(object, ...) {
if (object$P > 0)
{
cat("A univariate Normal probability distribution.\n")
cat("Expectation ", object$expectation, ", logged scale ", 
-0.5*log(object$P), ", standard deviation ", sqrt(1/object$P), ".\n", sep="")
}
else
cat("A flat univariate Normal probability distribution.\n")
}

simulate.normal <- function(object, nsim = 1, ...) {
   if (object$P>0) {
      rnorm(nsim, mean = object$expectation, sd = sqrt(1/object$P))
   }
   else {
      cat("ERROR: Cannot simulate from improper distribution.\n")
      rep(0, nsim)
   }
}

difference.normal <- function(object1, object2) {
if (class(object2)[1]=="normal") {
    normal(object2$expectation - object1$expectation, P = 1/(1/object1$P + 1/object2$P))
} else if (class(object2)[1]=="tdistribution") {
    cat("WARNING: The resulting distribution is an approximation.\n")
    tdistribution(object2$expectation - object1$expectation, 
    object2$degreesoffreedom*(1+object2$P/object1$P)^2, 
    P = 1/(1/object1$P + 1/object2$P))
} else
cat("ERROR: The difference function is not implemented for this pair of objects.\n")
}

linearpredict.normal <- function(object, X = 1 , PP = diag(length(X)), ...) {
   s <- length(X)
   if (s==1) PP <- matrix(PP, 1, 1)
   if (dim(PP)[1] != s | dim(PP)[2] != s) 
      cat("ERROR: P has wrong format.\n")
   else {
      expectation <- c(X*object$expectation, object$expectation)
      P   <- matrix(NA, s+1, s+1)
      P[1:s,1:s] <- PP
      P[1:s,s+1] <- -PP%*%X
      P[s+1,1:s] <- -PP%*%X
      P[s+1,s+1] <- X%*%PP%*%X + object$P
      result <- list(expectation=expectation, P=P)
      class(result) <- c("mnormal", "probabilitydistribution")
      result
   }
}

credibilityinterval.normal <- function(object, prob=0.95) {
c(qnorm((1-prob)/2, object$expectation, 1/sqrt(object$P)), 
  qnorm(1-(1-prob)/2, object$expectation, 1/sqrt(object$P)))
}

p.value.normal <- function(object, point=0) {
# For any distribution, this should return the probability outside the 
# credibility region with boundary containing the point. 
if (object$P==0) 
    1 
else
    pnorm(point, mean=object$expectation, sd = sqrt(1/object$P), 
          lower.tail = (object$expectation > point))*2
}

compose.normal <- function (object, type, ...) 
{
    if (type != "normal") 
        cat("ERROR: Not implemented for this type.\n")
    else {
        args <- list(...)
        mnormal(c(object$expectation, object$expectation), 
        matrix(c(object$P + args$P, -args$P, -args$P, args$P), 2, 2))
    }
}

##################################
## CLASS mnormal
##################################
# Data: 
# expectation: must be of length k>=1
# OPTIONAL: names: character vector of length k, giving the names of the dimensions
# P: must be symmetric matrix of size kxk (NEED NOT BE INVERTIBLE) 

mnormal <- function(expectation = c(0,0), P = diag(length(expectation))) {
   if (dim(P)[1]!=length(expectation))
      cat("ERROR in specification.\n")
   else if (dim(P)[2]!=length(expectation))
      cat("ERROR in specification.\n")
   else if (any(t(P)!=P)) 
      cat("ERROR: The precision matrix must be symmetric.\n")
   else if (min(eigen(P)$values)<0)
      cat("ERROR: The precision matrix must be non-negative definite.\n")
   else {
      result <- list(expectation=expectation, P=P)
      class(result) <- c("mnormal", "probabilitydistribution")
      result
   }
}

expectation.mnormal <- function(object) {
object$expectation
}

variance.mnormal <- function(object) {
   if (det(object$P)==0) cat("WARNING: The variance is infinite.\n")
   result <- ginv(object$P)
   0.5*(result + t(result))
}

precision.mnormal <- function(object) {
   object$P
}

marginal.mnormal <- function(object, v) { 
   # v gives the indices for the dimensions which should be kept. 
   k <- length(object$expectation)
   v <- as.integer(v)
   if (any(duplicated(v)))
      cat("ERROR: Duplicated indices not allowed.\n")
   else if (max(v)>k)
      cat("ERROR: Too large index.\n")
   else if (min(v)<1)
      cat("ERROR: Too small index.\n")
   else if (length(v)==k) {
      # permute the indices
      object$expectation <- object$expectation[v]
      object$P <- object$P[v,v]
      object
   }
   else {
      object$expectation <- object$expectation[v]
      object$P   <- object$P[v,v] - object$P[v,-v] %*% ginv(object$P[-v,-v]) %*% object$P[-v, v]
          object$P <- 0.5*(object$P + t(object$P))
      if (length(object$expectation)==1)
         class(object)[1] <- "normal"
      object
   }
}

conditional.mnormal <- function(object, v, val) {
# v gives the indices for the dimensions whose values should be fixed, and the 
# fixed values are given by val. 
   k <- length(object$expectation)
   v <- as.integer(v)
   if (any(duplicated(v)))
      cat("ERROR: Duplicated indices not allowed.\n")
   else if (max(v)>k)
      cat("ERROR: Too large index.\n")
   else if (min(v)<1)
      cat("ERROR: Too small index.\n")
   else if (length(v)==k)
      cat("ERROR: Cannot condition all variables.\n")
   else if (length(val)!=length(v))
      cat("ERROR: The length of the val vector must be equal to the number of indices.\n")
   else {
	object$expectation <- object$expectation[-v] - ginv(object$P[-v,-v])%*%
		object$P[-v,v]%*%(val-object$expectation[v])
	object$P <- object$P[-v,-v]
        if (length(object$expectation)==1)
           class(object)[1] <- "normal"
        object
   }
}

probabilitydensity.mnormal <- function(object, val, log=FALSE, normalize=TRUE) {
   result <- 0.5*log(det(object$P/(2*pi)))-0.5*(val-object$expectation)%*%object$P%*%(val-object$expectation)
   ifelse(log, result, exp(result))
}

print.mnormal <- function(x, ...) {
cat("A multivariate Normal probability distribution.\n")
cat("Expectation: ", x$expectation, "\n")
cat("Precision matrix:\n")
print(x$P) 
}

plot.mnormal <- function(x, onlybivariate=FALSE, ...) {
   if (onlybivariate & length(x$expectation)==2) {
      # To plot improper distributions reasonably: 
      Ptmp <- eigen(x$P)
      v <- ifelse(Ptmp$values, Ptmp$values, 1)
      Pnew <- Ptmp$vectors%*%diag(v)%*%t(Ptmp$vectors)
      st1 <- sqrt(Pnew[1,1]-Pnew[1,2]*Pnew[2,1]/Pnew[2,2])
      st2 <- sqrt(Pnew[2,2]-Pnew[2,1]*Pnew[1,2]/Pnew[1,1])
      startA <- x$expectation[1] - 3/st1
      stopA  <- x$expectation[1] + 3/st1
      startB <- x$expectation[2] - 3/st2
      stopB  <- x$expectation[2] + 3/st2
      logf <- function(X, ...) {
         -0.5*(X-x$expectation)%*%x$P%*%(X-x$expectation)
      }
      densitycontour(logf, c(startA, stopA, startB, stopB), data=0)
   } else {
      k <- length(x$expectation)
      old.par <- par(no.readonly = TRUE); on.exit(par(old.par))
      par(mfcol=c(k,k))
      for (i in 1:k) for (j in 1:k) 
         if (i==j) plot(marginal(x, i))
         else plot(marginal(x, c(i,j)), TRUE)
   }
}

summary.mnormal <- function(object, ...) {
cat("A multivariate Normal probability distribution.\n")
cat("Expectation ", object$expectation, "\n")
cat("Precision matrix:\n")
print(object$P) 
}

simulate.mnormal <- function(object, nsim = 1, ...) {
	if (det(object$P)==0) {
		cat("ERROR: Cannot simulate from improper distribution")
	} else {
		k <- length(object$expectation)
		B <- matrix(rnorm(nsim*k), nsim, k)%*%chol(solve(object$P))
		B + matrix(object$expectation, nsim, k, byrow=TRUE)
	}
}

linearpredict.mnormal <- function(object, X = rep(1,length(object$expectation)), PP = diag(length(X)/length(object$expectation)), ...) {
   if (is.matrix(X) && dim(X)[2] != length(object$expectation))
      cat("ERROR: X has wrong number of columns.\n")
   else if (!is.matrix(X) && length(X) != length(object$expectation))
      cat("ERROR: X has wrong length.\n")
   else if (dim(PP)[1] != length(X)/length(object$expectation) | dim(PP)[2] != length(X)/length(object$expectation))
      cat("ERROR: P as wrong format.\n")
   else {
      k <- length(object$expectation)
      s <- length(X)/k
      X <- matrix(X, s, k)
      expectation <- c(X%*%object$expectation, object$expectation)
      P   <- matrix(NA, s+k, s+k)
      P[1:s,1:s] <- PP
      P[1:s,s+(1:k)] <- -PP%*%X
      P[s+(1:k),1:s] <- -t(PP%*%X)
      P[s+(1:k),s+(1:k)] <- t(X)%*%PP%*%X + object$P
      result <- list(expectation=expectation, P=P)
      class(result) <- c("mnormal", "probabilitydistribution")
      result
   }
}

p.value.mnormal <- function(object, point=0) {
# For any distribution, this should return the probability outside the 
# credibility region with boundary containing the origin. 
    if (sum(point*point)==0) point <- rep(0, length(object$expectation))
    pchisq((object$expectation-point)%*%object$P%*%(object$expectation-point) ,
		length(object$expectation), lower.tail=FALSE)
}


##################################
## CLASS gammadistribution 
##################################
# Data: 
# alpha >= 0
# beta >= 0

gammadistribution <- function(alpha = 1, beta = 1) {
if (alpha<0) 
cat("ERROR in specification\n")
else if (beta<0)
cat("ERROR in specification\n")
else {
	result <- list(alpha=alpha, beta=beta)
	class(result) <- c("gammadistribution", "probabilitydistribution")
	result
}
}

expectation.gammadistribution <- function(object) {
object$alpha/object$beta
}

variance.gammadistribution <- function(object) {
object$alpha/object$beta^2
}

precision.gammadistribution <- function(object) {
object$beta^2/object$alpha
}

probabilitydensity.gammadistribution <- function(object, val, log=FALSE, normalize=TRUE) {
   dgamma(val, object$alpha, object$beta, log=log)
}

cdf.gammadistribution <- function(object, val) {
   pgamma(val, object$alpha, object$beta)
}

invcdf.gammadistribution <- function(object, val) {
   if (val<=0 | val>=1) 
      cat("ERROR: Illegal input for the invcdf function.\n")
   else 
      qgamma(val, object$alpha, object$beta)
}

print.gammadistribution <- function(x, ...) {
cat("A Gamma probability distribution.\n")
cat("Parameters: ", x$alpha, x$beta, "\n")
cat("Expectation: ", x$alpha/x$beta, "\n")
cat("Variance: ", x$alpha/x$beta^2, "\n")
}

plot.gammadistribution <- function(x, ...) {
   start <- qgamma(0.000001, x$alpha, x$beta)
   stop  <- qgamma(0.999, x$alpha, x$beta)
   xx    <- seq(start, stop, length=1001)
   yy    <- dgamma(xx, x$alpha, x$beta)
   plot(xx, yy, xlab = "", ylab="", type="l")
}

summary.gammadistribution <- function(object, ...) {
cat("A Gamma probability distribution.\n")
cat("Parameters: ", object$alpha, object$beta, "\n")
cat("Expectation: ", object$alpha/object$beta, "\n")
cat("Variance: ", object$alpha/object$beta^2, "\n")
}

simulate.gammadistribution <- function(object, nsim=1, ...) {
   rgamma(nsim, object$alpha, object$beta)
}

p.value.gammadistribution <- function(object, point=0) {
# For any distribution, this should return the probability outside the 
# credibility region with boundary containing the point. 
    res <- pgamma(point, object$alpha, object$beta)
    if (res>0.5) res <- 1 -res
    2*res
}

credibilityinterval.gammadistribution <- function(object, prob=0.95) {
c(qgamma((1-prob)/2, object$alpha, object$beta), 
  qgamma(1-(1-prob)/2, object$alpha, object$beta))
}

##################################
## CLASS expgamma 
##################################
# Data: 
# alpha >= 0
# beta >= 0
# gamma

expgamma <- function(alpha = 1, beta = 1, gamma = -2) {
if (alpha<0) 
cat("ERROR in specification\n")
else if (beta<0)
cat("ERROR in specification\n")
else {
	result <- list(alpha=alpha, beta=beta, gamma=gamma)
	class(result) <- c("expgamma", "probabilitydistribution")
	result
}
}

expectation.expgamma <- function(object) {
(digamma(object$alpha) - log(object$beta))/object$gamma
}

variance.expgamma <- function(object) {
trigamma(object$alpha)/object$gamma^2
}

precision.expgamma <- function(object) {
object$gamma^2/trigamma(object$alpha)
}

probabilitydensity.expgamma <- function(object, val, log=FALSE, normalize=TRUE) {
   x <- exp(object$gamma*val)
   dgamma(x, object$alpha, object$beta)*x*abs(object$gamma)
}

cdf.expgamma <- function(object, val) {
   x <- exp(object$gamma*val)
   pgamma(x, object$alpha, object$beta, lower.tail=(object$gamma>0))
}

invcdf.expgamma <- function(object, val) {
   if (val<=0 | val>=1) 
      cat("ERROR: Illegal input for the invcdf function.\n")
   else 
      log(qgamma(val, object$alpha, object$beta, lower.tail=(object$gamma>0)))/object$gamma
}

print.expgamma <- function(x, ...) {
cat("An ExpGamma probability distribution.\n")
cat("Parameters: ", x$alpha, x$beta, x$gamma, "\n")
}

plot.expgamma <- function(x, ...) {
   start <- invcdf(x, 0.001)
   stop  <- invcdf(x, 0.999)
   xx    <- seq(start, stop, length=1001)
   yy    <- probabilitydensity(x, xx)
   plot(xx, yy, xlab = "", ylab="", type="l")
}

summary.expgamma <- function(object, ...) {
cat("An ExpGamma probability distribution.\n")
cat("Parameters: ", object$alpha, object$beta, object$gamma, "\n")
}

simulate.expgamma <- function(object, nsim=1, ...) {
   log(rgamma(nsim, object$alpha, object$beta))/object$gamma
}

p.value.expgamma <- function(object, point=0) {
# For any distribution, this should return the probability outside the 
# credibility region with boundary containing the point. 
    res <- pgamma(exp(object$gamma*point), object$alpha, object$beta)
    if (res>0.5) res <- 1 -res
    2*res
}

credibilityinterval.expgamma <- function(object, prob=0.95) {
   c(invcdf(object, (1-prob)/2), invcdf(object, 1-(1-prob)/2))
}

##################################
## CLASS normalgamma 
##################################
# Data: 
# OPTIONAL: name: giving the names of the dimension
# mu 
# P >= 0: Although externally, we use the parameter kappa = -0.5*log(P) 
# alpha >= 0
# beta >= 0

normalgamma <- function(mu, kappa , alpha, beta) {
if (missing(mu) & missing(kappa) & missing(alpha) & missing(beta)) {
   # The standard prior: 
   result <- list(mu = 0, P = 0, alpha=-0.5, beta=0)
   class(result) <- c("normalgamma", "probabilitydistribution")
   result 
}
else if (length(mu)!=1)
   cat("ERROR: A single number only for the first parameter.\n")
else {
   result <- list(mu = mu, P = exp(-2*kappa), alpha=alpha, beta=beta)
   class(result) <- c("normalgamma", "probabilitydistribution")
   result
}
}

expectation.normalgamma <- function(object) {
c(expectation(marginal(object, 1)), expectation(marginal(object, 2)))
}

variance.normalgamma <- function(object) {
diag(variance(marginal(object, 1)), variance(marginal(object,2)))
}

precision.normalgamma <- function(object) {
diag(precision(marginal(object, 1)), precision(marginal(object, 2)))
}

marginal.normalgamma <- function(object, v) {
# v gives the index for the dimension which should be kept. 
   v <- as.integer(v)
   if (length(v)>1) 
      cat("ERROR: Too long vector of indices for the marginal.\n")
   else if (v>2)
      cat("ERROR: Too large index.\n")
   else if (v<1)
      cat("ERROR: Too small index.\n")
   else {
      if (v==1) 
         tdistribution(object$mu, 2*object$alpha, P = object$alpha/object$beta*object$P)
      else 
         gammadistribution(object$alpha, object$beta)
   }
}

conditional.normalgamma <- function(object, v, val) {
# v gives the indices for the dimensions whose values should be fixed, and the 
# fixed values are given by val. 
   v <- as.integer(v)
   if (length(v)>1)
      cat("ERROR: Too long vector of indices for the conditional.\n")
   else if (v>2)
      cat("ERROR: Too large index.\n")
   else if (v<1)
      cat("ERROR: Too small index.\n")
   else if (length(val)>1)
      cat("ERROR: The number of values must be equal to the number of indices.\n")
   else {
      if (v==1) 
         gammadistribution(object$alpha+0.5, object$beta + object$P/2*(val-object$mu)^2)
      else 
         normal(object$mu, P=object$P*val)
   }
}

probabilitydensity.normalgamma <- function(object, val, log=FALSE, normalize=TRUE) {
   result <- dgamma(val[2], object$alpha, object$beta) * dnorm(val[1], object$mu, 1/sqrt(val[2]*object$P))
   ifelse(log, log(result), result)
}

print.normalgamma <- function(x, ...) {
cat("A bivariate Normal-Gamma probability distribution.\n")
cat("mu:    ", x$mu, "\n")
cat("kappa: ", -0.5*log(x$P),  "\n")
cat("alpha: ", x$alpha, "\n")
cat("beta:  ", x$beta, "\n")
}

plot.normalgamma <- function(x, onlybivariate=FALSE, switchaxes=FALSE, ...) {
# The onlybivariate parameter can make the output consist of only the bivariate plot
# The switchaxes parameter can make the axes switch, if the plot is only bivariate. 
   if (onlybivariate) {
      A <- marginal(x, 1)
      B <- marginal(x, 2)
      startA <- qt(0.001, A$degreesoffreedom)/sqrt(A$P) + A$expectation
      stopA  <- qt(0.999, A$degreesoffreedom)/sqrt(A$P) + A$expectation
      startB <- qgamma(0.000001, B$alpha, B$beta)
      stopB  <- qgamma(0.999, B$alpha, B$beta)
      logf <- function(XX, ...) {
         (x$alpha - 0.5)*log(XX[2]) - XX[2]*(x$beta + x$P/2*(XX[1]-x$mu)^2)
      }
      logft <- function(XX, ...) {
         (x$alpha - 0.5)*log(XX[1]) - XX[1]*(x$beta + x$P/2*(XX[2]-x$mu)^2)
      }
      if (switchaxes) 
         densitycontour(logf, c(startA, stopA, startB, stopB), data=0)
      else
         densitycontour(logft, c(startB, stopB, startA, stopA), data=0)   
   } else if (x$alpha==-0.5 & x$beta==0.0 & x$mu==0.0 & x$P==0.0) {
      cat("A flat distribution.\n")
   } else {
      old.par <- par(no.readonly = TRUE); on.exit(par(old.par))
      par(mfcol=c(2,2))
      plot(marginal(x, 1))
      plot(x, TRUE, TRUE)
      plot(x, TRUE, FALSE)
      plot(marginal(x, 2))
   }
}

summary.normalgamma <- function(object, ...) {
print(object)
}

simulate.normalgamma <- function(object, nsim=1, ...) {
   tau <- rgamma(nsim, object$alpha, object$beta)
   X   <- rnorm(nsim)
   X   <- X/sqrt(tau*object$P) + object$mu
   cbind(X, tau)
}

linearpredict.normalgamma <- function(object, X = 1 , P = diag(length(as.vector(X))), ...) {
   X <- as.vector(X)
   s <- length(X)
   if (dim(P)[1] != s | dim(P)[2] != s) 
      cat("ERROR: P has wrong format.\n")
   else {
      expectation <- c(X*object$mu, object$mu)
      precision   <- matrix(NA, s+1, s+1)
      precision[1:s,1:s] <- P
      precision[1:s,s+1] <- -P%*%X
      precision[s+1,1:s] <- -P%*%X
      precision[s+1,s+1] <- X%*%P%*%X + object$P
      result <- list(mu=expectation, P=precision, alpha=object$alpha, beta=object$beta)
      class(result) <- c("mnormalgamma", "probabilitydistribution")
      result
   }
}

prediction.normalgamma <- function(object) {
# This should give the distribution when making a NEW observation with 
# the theta of this distribution as expectation, and the tau of this 
# distribution as precision
   tdistribution(object$mu, 2*object$alpha, P = object$alpha/object$beta*object$P/(object$P + 1))
}

##################################
## CLASS normalexpgamma 
##################################
# Data: 
# OPTIONAL: name: giving the names of the dimension
# mu 
# P >= 0: Although externally, we use the parameter kappa = -0.5*log(P) 
# alpha >= 0
# beta >= 0

normalexpgamma <- function(mu, kappa, alpha, beta) {
if (missing(mu) & missing(kappa) & missing(alpha) & missing(beta)) {
   # The standard prior: 
   result <- list(mu = 0, P = 0, alpha=-0.5, beta=0)
   class(result) <- c("normalexpgamma", "probabilitydistribution")
   result 
}
else if (length(mu)!=1)
   cat("ERROR: A single number only for the first parameter.\n")
else {
   result <- list(mu = mu, P = exp(-2*kappa), alpha=alpha, beta=beta)
   class(result) <- c("normalexpgamma", "probabilitydistribution")
   result
}
}

expectation.normalexpgamma <- function(object) {
c(expectation(marginal(object, 1)), expectation(marginal(object, 2)))
}

variance.normalexpgamma <- function(object) {
diag(c(variance(marginal(object, 1)), variance(marginal(object, 2))))
}

precision.normalexpgamma <- function(object) {
diag(c(precision(marginal(object, 1)), precision(marginal(object, 2))))
}

marginal.normalexpgamma <- function(object, v) {
# v gives the index for the dimension which should be kept. 
   v <- as.integer(v)
   if (length(v)>1) 
      cat("ERROR: Too long vector of indices for the marginal.\n")
   else if (v>2)
      cat("ERROR: Too large index.\n")
   else if (v<1)
      cat("ERROR: Too small index.\n")
   else {
      if (v==1) 
         tdistribution(object$mu, 2*object$alpha, P = object$alpha/object$beta*object$P)
      else 
         expgamma(object$alpha, object$beta, -2)
   }
}

conditional.normalexpgamma <- function(object, v, val) {
# v gives the indices for the dimensions whose values should be fixed, and the 
# fixed values are given by val. 
   v <- as.integer(v)
   if (length(v)>1)
      cat("ERROR: Too long vector of indices for the conditional.\n")
   else if (v>2)
      cat("ERROR: Too large index.\n")
   else if (v<1)
      cat("ERROR: Too small index.\n")
   else if (length(val)>1)
      cat("ERROR: The number of values must be equal to the number of indices.\n")
   else {
      if (v==1) 
         expgamma(object$alpha+0.5, object$beta + object$P/2*(val-object$mu)^2, -2)
      else 
         normal(object$mu, P=object$P*exp(-2*val))
   }
}

probabilitydensity.normalexpgamma <- function(object, val, log=FALSE, normalize=TRUE) {
   x <- exp(-2*val[2])
   result <- dgamma(x, object$alpha, object$beta) * dnorm(val[1], object$mu, 1/sqrt(x*object$P))*x*2
   ifelse(log, log(result), result)
}

print.normalexpgamma <- function(x, ...) {
cat("A bivariate Normal-ExpGamma probability distribution.\n")
cat("mu:    ", x$mu, "\n")
cat("kappa: ", -0.5*log(x$P),  "\n")
cat("alpha: ", x$alpha, "\n")
cat("beta:  ", x$beta, "\n")
}

plot.normalexpgamma <- function(x, onlybivariate=FALSE, switchaxes=FALSE, ...) {
# The onlybivariate parameter can make the output consist of only the bivariate plot
# The switchaxes parameter can make the axes switch, if the plot is only bivariate. 
   if (x$mu==0 & x$P==0 & x$alpha==-0.5 & x$beta==0) 
      cat("The standard flat improper prior.\n")
   else if (onlybivariate) {
      A <- marginal(x, 1)
      B <- marginal(x, 2)
      startA <- invcdf(A, 0.001)
      stopA  <- invcdf(A, 0.999)
      startB <- invcdf(B, 0.001)
      stopB  <- invcdf(B, 0.999)
      logf <- function(XX, ...) {
         - (2*x$alpha + 1)*XX[2] - exp(-2*XX[2])*(x$beta + x$P/2*(XX[1]-x$mu)^2)
      }
      logft <- function(XX, ...) {
         - (2*x$alpha + 1)*XX[1] - exp(-2*XX[1])*(x$beta + x$P/2*(XX[2]-x$mu)^2)
      }
      if (switchaxes) 
         densitycontour(logf, c(startA, stopA, startB, stopB), data=0)
      else
         densitycontour(logft, c(startB, stopB, startA, stopA), data=0)   
   } else {
      old.par <- par(no.readonly = TRUE); on.exit(par(old.par))
      par(mfcol=c(2,2))
      plot(marginal(x, 1))
      plot(x, TRUE, TRUE)
      plot(x, TRUE, FALSE)
      plot(marginal(x, 2))
   }
}

summary.normalexpgamma <- function(object, ...) {
print(object)
}

simulate.normalexpgamma <- function(object, nsim=1, ...) {
   tau <- rgamma(nsim, object$alpha, object$beta)
   X   <- rnorm(nsim)
   X   <- X/sqrt(tau*object$P) + object$mu
   cbind(X, -0.5*log(tau))
}

linearpredict.normalexpgamma <- function(object, X = 1 , P = diag(length(as.vector(X))), ...) {
   X <- as.vector(X)
   s <- length(X)
   if (dim(P)[1] != s | dim(P)[2] != s) 
      cat("ERROR: P has wrong format.\n")
   else {
      expectation <- c(X*object$mu, object$mu)
      precision   <- matrix(NA, s+1, s+1)
      precision[1:s,1:s] <- P
      precision[1:s,s+1] <- -P%*%X
      precision[s+1,1:s] <- -P%*%X
      precision[s+1,s+1] <- X%*%P%*%X + object$P
      result <- list(mu=expectation, P=precision, alpha=object$alpha, beta=object$beta)
      class(result) <- c("mnormalexpgamma", "probabilitydistribution")
      result
   }
}

prediction.normalexpgamma <- function(object) {
# This should give the distribution when making a NEW observation with 
# the theta of this distribution as expectation, and the tau of this 
# distribution as precision
     tdistribution(object$mu, 2*object$alpha, P = object$alpha/object$beta*object$P/(object$P + 1))
}


##################################
## CLASS mnormalgamma 
##################################
# Data: 
# OPTIONAL: names: character vector of length k, giving the names of the dimensions
# mu 
# P non-negative symmetric matrix
# alpha >= 0
# beta >= 0


## CONSTRUCTORS: 
################
mnormalgamma <- function(mu=c(0,0), P , alpha, beta) {
if (missing(P) & missing(alpha) & missing(beta)) {
   #Use the standard prior: 
   k <- length(mu)
   result <- list(mu = mu, P = matrix(0, k, k), alpha=-k/2, beta=0)
   class(result) <- c("mnormalgamma", "probabilitydistribution")
   result
}
else {
   if (missing(P)) P <- diag(length(mu))
   if (any(eigen(P)$values<0)) 
      cat("ERROR: Illegal value for parameter P.\n")
   else {    
      result <- list(mu = mu, P = P, alpha=alpha, beta=beta)
      class(result) <- c("mnormalgamma", "probabilitydistribution")
      result
   }
}
}

## Generic functions: 
######################

expectation.mnormalgamma <- function(object) {
   c(object$mu, object$alpha/object$beta)
}

variance.mnormalgamma <- function(object, ...) {
   k <- length(object$mu)
   result <- matrix(0, k+1, k+1)
   result[1:k,1:k] <- object$beta/(object$alpha-1)*ginv(object$P)
   result[k+1,k+1] <- object$alpha/object$beta^2
   result <- 0.5*(result + t(result))
   result
}

precision.mnormalgamma <- function(object, ...) {
   k <- length(object$mu)
   result <- matrix(0, k+1, k+1)
   result[1:k,1:k] <- (object$alpha-1)/object$beta*object$P
   result[k+1,k+1] <- object$beta^2/object$alpha
   result
}

marginal.mnormalgamma <- function(object, v) {
   # v gives the indices for the dimensions which should be kept. 
   k <- length(object$mu)
   v <- as.integer(v)
   if (any(duplicated(v)))
      cat("ERROR: Duplicated indices not allowed.\n")
   else if (max(v)>k+1)
      cat("ERROR: Too large index.\n")
   else if (min(v)<1)
      cat("ERROR: Too small index.\n")
   else {
      KeepGamma <- (max(v)==k+1)
      w <- v[v<=k]
      if (length(w)==0) 
         gammadistribution(object$alpha, object$beta)
      else {
         object$mu <- object$mu[w]
         if (length(w)==k)
            object$P <- object$P[w,w]
         else {
            object$P <- object$P[w,w] - object$P[w,-w]%*%ginv(object$P[-w,-w])%*%object$P[-w,w]
            object$P <- 0.5*(object$P + t(object$P))
         }
         if (KeepGamma) {
            if (length(object$mu)>1)
               object
            else 
               normalgamma(object$mu, -0.5*log(object$P), object$alpha, object$beta)
         } else {
            if (length(object$mu)>1)
    	       mtdistribution(object$mu, 2*object$alpha, P = 
                  rep(object$alpha/object$beta, prod(dim(object$P)))*object$P) 
            else
      	       tdistribution(object$mu, 2*object$alpha, P = object$alpha/object$beta*object$P) 
         }
      }
   }
}

conditional.mnormalgamma <- function(object, v, val) {
# v gives the indices for the dimensions whose values should be fixed, and the 
# fixed values are given by val.
   if (missing(v)) 
      cat("ERROR: Argument giving indices for conditioning is missing.\n")
   else if (missing(val))
      cat("ERROR: Argument giving values for conditioning is missing.\n")
   else { 
   k <- length(object$mu)
   v <- as.integer(v)
   if (length(v)<1) 
      cat("ERROR: Cannot condition on empty index vector.\n")
   else if (any(duplicated(v)))
      cat("ERROR: Duplicated indices not allowed.\n")
   else if (max(v)>k+1)
      cat("ERROR: Too large index.\n")
   else if (min(v)<1)
      cat("ERROR: Too small index.\n")
   else if (length(v)==k+1)
      cat("ERROR: Cannot condition on all indices simultaneously.\n")
   else if (length(val)!=length(v))
      cat("ERROR: The number of values must be equal to the number of indices.\n")
   else if (max(v)==k+1) {
      obj <- mnormal(object$mu, P=rep(val[v==k+1],prod(dim(object$P)))*object$P)
      if (length(v)>1) 
         conditional(obj, v[v<=k], val[v<=k])
      else 
         obj
   } else if (length(v)==k) {
      gammadistribution(object$alpha + k/2, object$beta + 
         0.5*(val-object$mu)%*%object$P%*%(val-object$mu))
   } else {
      newexp <- as.vector(object$mu[-v] - ginv(object$P[-v,-v])%*%object$P[-v,v]%*%(val - object$mu[v]))
      newbeta <- object$beta + 0.5*(val - object$mu[v])%*%
         (object$P[v,v]-object$P[v,-v]%*%ginv(object$P[-v,-v])%*%object$P[-v,v])%*%
         (val - object$mu[v])
      if (length(newexp)>1) {
         object$mu <- newexp
         object$P  <- object$P[-v,-v]
	 object$alpha <- object$alpha + length(v)/2
         object$beta <- newbeta
         object
      } else 
         normalgamma(newexp, -0.5*log(object$P[-v,-v]), object$alpha + length(v)/2, newbeta)
   }      
   }
}

probabilitydensity.mnormalgamma <- function(object, val, log=FALSE, normalize=TRUE) {
   k <- length(object$mu)
   result <- dgamma(val[k+1], object$alpha, object$beta, log=TRUE)  
   0.5*log(det(object$P*val[k+1]/(2*pi)))-0.5*(val[1:k]-object$mu)%*%object$P%*%(val[1:k]-object$mu)
   ifelse(log, result, exp(result))
}

print.mnormalgamma <- function(x, ...) {
cat("A Multivariate Normal-Gamma probability distribution.\n")
cat("mu:    ", x$mu, "\n")
cat("P:     ", x$P,  "\n")
cat("alpha: ", x$alpha, "\n")
cat("beta:  ", x$beta, "\n")
}

plot.mnormalgamma <- function(x, ...) {
   if (max(abs(x$mu))==0 & max(abs(x$P))==0 & x$alpha==-length(x$mu)/2 & x$beta==0) {
      cat("A flat prior.\n")
   } else {
   k <- length(x$mu) + 1
   old.par <- par(no.readonly = TRUE); on.exit(par(old.par))
   par(mfcol=c(k,k))
   for (i in 1:k) for (j in 1:k) {
      if (i==j)
         plot(marginal(x, i))
      else if (j==k) 
         plot(marginal(x, c(i,j)), TRUE, TRUE)
      else
         plot(marginal(x, c(i,j)), TRUE) 
   }
}
}

summary.mnormalgamma <- function(object, ...) {
print(object)
}

simulate.mnormalgamma <- function(object, nsim=1, ...) {
   tau <- rgamma(nsim, object$alpha, object$beta)
   k   <- length(object$mu)
   X   <- matrix(rnorm(nsim * k), nsim, k) %*% chol(solve(object$P))
   X   <- X * matrix(rep(1/sqrt(tau),k), nsim, k)
   X   <- X + matrix(object$mu, nsim, k, byrow=TRUE)
   cbind(X, tau)
}

linearpredict.mnormalgamma <- function(object, X = rep(1,length(object$mu)), P = diag(length(X)/length(object$mu)), ...) {
   if (is.matrix(X) && dim(X)[2] != length(object$mu))
      cat("ERROR: X has wrong number of columns.\n")
   else if (!is.matrix(X) && length(X) != length(object$mu))
      cat("ERROR: X has wrong length.\n")
   else if (dim(P)[1] != length(X)/length(object$mu) | dim(P)[2] != length(X)/length(object$mu))
      cat("ERROR: P as wrong format.\n")
   else {
      k <- length(object$mu)
      s <- length(X)/k
      X <- matrix(X, s, k)
      expectation <- c(X%*%object$mu, object$mu)
      precision   <- matrix(NA, s+k, s+k)
      precision[1:s,1:s] <- P
      precision[1:s,s+(1:k)] <- -P%*%X
      precision[s+(1:k),1:s] <- -t(P%*%X)
      precision[s+(1:k),s+(1:k)] <- t(X)%*%P%*%X + object$P
      result <- list(mu=expectation, P=precision, alpha=object$alpha, beta=object$beta)
      class(result) <- c("mnormalgamma", "probabilitydistribution")
      result
   }
}

contrast.mnormalgamma <- function(object, v) {
   if (length(object$mu)!=length(v))
      cat("ERROR: Wrong length of the input vector.\n")
   else 
      normalgamma(sum(v*object$mu), 0.5*log(as.vector(v%*%ginv(object$P)%*%v)), object$alpha, object$beta)
}

##################################
## CLASS mnormalexpgamma 
##################################
# Data: 
# OPTIONAL: names: character vector of length k, giving the names of the dimensions
# mu 
# P non-negative symmetric matrix
# alpha >= 0
# beta >= 0


## CONSTRUCTORS: 
################
mnormalexpgamma <- function(mu=c(0,0), P , alpha, beta) {
if (missing(P) & missing(alpha) & missing(beta)) {
   #Use the standard prior: 
   k <- length(mu)
   result <- list(mu = mu, P = matrix(0, k, k), alpha=-k/2, beta=0)
   class(result) <- c("mnormalexpgamma", "probabilitydistribution")
   result
}
else { 
   if (missing(P)) P <- diag(length(mu))   
   result <- list(mu = mu, P = P, alpha=alpha, beta=beta)
   class(result) <- c("mnormalexpgamma", "probabilitydistribution")
   result
}
}

## Generic functions: 
######################

expectation.mnormalexpgamma <- function(object) {
k <- length(object$mu)
c(object$mu, expectation(marginal(object, k+1)))
}

variance.mnormalexpgamma <- function(object, ...) {
k <- length(object$mu)
result <- matrix(0, k+1, k+1)
result[1:k,1:k] <- variance(marginal(object, 1:k))
result[k+1,k+1] <- variance(marginal(object, k+1))
result
}

precision.mnormalexpgamma <- function(object, ...) {
k <- length(object$mu)
result <- matrix(0, k+1, k+1)
result[1:k,1:k] <- precision(marginal(object, 1:k))
result[k+1,k+1] <- precision(marginal(object, k+1))
result
}

marginal.mnormalexpgamma <- function(object, v) {
# v gives the indices for the dimensions which should be kept. 
   k <- length(object$mu)
   v <- as.integer(v)
   if (any(duplicated(v)))
      cat("ERROR: Duplicated indices not allowed.\n")
   else if (max(v)>k+1)
      cat("ERROR: Too large index.\n")
   else if (min(v)<1)
      cat("ERROR: Too small index.\n")
   else {
      KeepGamma <- (max(v)==k+1)
      w <- v[v<=k]
      if (length(w)==0) 
         expgamma(object$alpha, object$beta)
      else {
        object$mu <- object$mu[w]
         if (length(w)==k)
            object$P <- object$P[w,w]
         else {
            object$P <- object$P[w,w] - object$P[w,-w]%*%ginv(object$P[-w,-w])%*%object$P[-w,w]
            object$P <- 0.5*(object$P + t(object$P))
         }
         if (KeepGamma) {
            if (length(object$mu)>1)
               object
            else 
               normalexpgamma(object$mu, -0.5*log(object$P), object$alpha, object$beta)
         } else {
            if (length(object$mu)>1)
    	       mtdistribution(object$mu, 2*object$alpha, P = 
                  rep(object$alpha/object$beta, prod(dim(object$P)))*object$P) 
            else
      	       tdistribution(object$mu, 2*object$alpha, P = object$alpha/object$beta*object$P) 
         }
      }
   }
}


conditional.mnormalexpgamma <- function(object, v, val) {
# v gives the indices for the dimensions whose values should be fixed, and the 
# fixed values are given by val.
   if (missing(v)) 
      cat("ERROR: Argument giving indices for conditioning is missing.\n")
   else if (missing(val))
      cat("ERROR: Argument giving values for conditioning is missing.\n")
   else { 
   k <- length(object$mu)
   v <- as.integer(v)
   if (length(v)<1) 
      cat("ERROR: Cannot condition on empty index vector.\n")
   else if (any(duplicated(v)))
      cat("ERROR: Duplicated indices not allowed.\n")
   else if (max(v)>k+1)
      cat("ERROR: Too large index.\n")
   else if (min(v)<1)
      cat("ERROR: Too small index.\n")
   else if (length(v)==k+1)
      cat("ERROR: Cannot condition on all indices simultaneously.\n")
   else if (length(val)!=length(v))
      cat("ERROR: The number of values must be equal to the number of indices.\n")
   else if (max(v)==k+1) {
      obj <- mnormal(object$mu, P=rep(val[v==k+1],prod(dim(object$P)))*object$P)
      if (length(v)>1) 
         conditional(obj, v[v<=k], val[v<=k])
      else 
         obj
   } else if (length(v)==k) {
      expgamma(object$alpha + k/2, object$beta + 
         0.5*(val-object$mu)%*%object$P%*%(val-object$mu), -2)
   } else {
      newexp <- as.vector(object$mu[-v] - ginv(object$P[-v,-v])%*%object$P[-v,v]%*%(val - object$mu[v]))
      newbeta <- object$beta + 0.5*(val - object$mu[v])%*%
         (object$P[v,v]-object$P[v,-v]%*%ginv(object$P[-v,-v])%*%object$P[-v,v])%*%
         (val - object$mu[v])
      if (length(newexp)>1) {
         object$mu <- newexp
         object$P  <- object$P[-v,-v]
	 object$alpha <- object$alpha + length(v)/2
         object$beta <- newbeta
         object
      } else 
         normalexpgamma(newexp, -0.5*log(object$P[-v,-v]), object$alpha + length(v)/2, newbeta)
   }      
   }
}

probabilitydensity.mnormalexpgamma <- function(object, val, log=FALSE, normalize=TRUE) {
   k <- length(object$mu)
   result <- dgamma(exp(-2*val[k+1]), object$alpha, object$beta, log=TRUE)  
   0.5*log(det(object$P*val[k+1]/(2*pi)))-0.5*(val[1:k]-object$mu)%*%object$P%*%(val[1:k]-object$mu)
   ifelse(log, result, exp(result))
}

print.mnormalexpgamma <- function(x, ...) {
cat("A Multivariate Normal-ExpGamma probability distribution.\n")
cat("mu:    ", x$mu, "\n")
cat("P:     ", x$P,  "\n")
cat("alpha: ", x$alpha, "\n")
cat("beta:  ", x$beta, "\n")
}

plot.mnormalexpgamma <- function(x, ...) {
   if (max(abs(x$mu))==0 & max(abs(x$P))==0 & x$alpha==-length(x$mu)/2 & x$beta==0) {
      cat("A flat prior.\n")
   } else {
   k <- length(x$mu) + 1
   old.par <- par(no.readonly = TRUE); on.exit(par(old.par))
   par(mfcol=c(k,k))
   for (i in 1:k) for (j in 1:k) {
      if (i==j)
         plot(marginal(x, i))
      else if (j==k) 
         plot(marginal(x, c(i,j)), TRUE, TRUE)
      else
         plot(marginal(x, c(i,j)), TRUE) 
   }
}
}

summary.mnormalexpgamma <- function(object, ...) {
print(object)
}

simulate.mnormalexpgamma <- function(object, nsim=1, ...) {
   tau <- rgamma(nsim, object$alpha, object$beta)
   k   <- length(object$mu)
   X   <- matrix(rnorm(nsim * k), nsim, k) %*% chol(ginv(object$P))
   X   <- X * matrix(rep(1/sqrt(tau),k), nsim, k)
   X   <- X + matrix(object$mu, nsim, k, byrow=TRUE)
   cbind(X, -0.5*log(tau))
}

linearpredict.mnormalexpgamma <- function(object, X = rep(1,length(object$mu)), P = diag(length(X)/length(object$mu)), ...) {
   if (is.matrix(X) && dim(X)[2] != length(object$mu))
      cat("ERROR: X has wrong number of columns.\n")
   else if (!is.matrix(X) && length(X) != length(object$mu))
      cat("ERROR: X has wrong length.\n") 
   else if (dim(P)[1] != length(X)/length(object$mu) | dim(P)[2] != length(X)/length(object$mu))
      cat("ERROR: P as wrong format.\n")
   else {
      k <- length(object$mu)
      s <- length(X)/k
      X <- matrix(X, s, k)
      expectation <- c(X%*%object$mu, object$mu)
      precision   <- matrix(NA, s+k, s+k)
      precision[1:s,1:s] <- P
      precision[1:s,s+(1:k)] <- -P%*%X
      precision[s+(1:k),1:s] <- -t(P%*%X)
      precision[s+(1:k),s+(1:k)] <- t(X)%*%P%*%X + object$P
      result <- list(mu=expectation, P=precision, alpha=object$alpha, beta=object$beta)
      class(result) <- c("mnormalexpgamma", "probabilitydistribution")
      result
   }
}

contrast.mnormalexpgamma <- function(object, v) {
   if (length(object$mu)!=length(v))
      cat("ERROR: Wrong length of the input vector.\n")
   else 
      normalexpgamma(sum(v*object$mu), 0.5*log(as.vector(v%*%ginv(object$P)%*%v)), object$alpha, object$beta)
}

##################################
## CLASS tdistribution 
##################################
# Data: 
# expectation
# degreesoffreedom
# P

tdistribution <- function(expectation=0, degreesoffreedom = 1e20, lambda, P = 1) {
if (!missing(lambda) & !missing(P))
   cat("ERROR: Cannot specify both logged scale and P at the same time.")
else if (degreesoffreedom < 1) 
  cat("ERROR: The degrees of freedom must be at least 1.\n") 
else if (P < 0 )
  cat("ERROR: The P parameter cannot be negative.\n")
else {
   if (!missing(lambda)) P <- exp(-2*lambda)
result <- list(expectation=as.vector(expectation), degreesoffreedom=degreesoffreedom, P=as.vector(P))
class(result) <- c("tdistribution", "probabilitydistribution")
result
}
}

expectation.tdistribution <- function(object) {
object$expectation
}

variance.tdistribution <- function(object) {
if (object$degreesoffreedom>2) 
   object$degreesoffreedom/(object$degreesoffreedom-2)/object$P
else
   1/0
} 

precision.tdistribution <- function(object) {
object$precision
if (object$degreesoffreedom>2) 
   1/(object$degreesoffreedom/(object$degreesoffreedom-2)/object$P)
else
   1/0
}

probabilitydensity.tdistribution <- function(object, val, log=FALSE, normalize=TRUE) {
   dt((val - object$expectation)*sqrt(object$P), object$degreesoffreedom, log=log)
}

cdf.tdistribution <- function(object, val) {
   pt((val-object$expectation)*sqrt(object$P), object$degreesoffreedom)
}

invcdf.tdistribution <- function(object, val) {
   if (val<=0 | val>=1) 
      cat("ERROR: Illegal input for the invcdf function.\n")
   else 
      qt(val, object$degreesoffreedom)/sqrt(object$P) + object$expectation
}

print.tdistribution <- function(x, ...) {
cat("A Student t probability distribution.\n")
cat("Expectation ", x$expectation, ", degrees of freedom ", x$degreesoffreedom, 
    ", logged scale ", -0.5*log(x$P), "\n")
}

plot.tdistribution <- function(x, ...) {
   start <- qt(0.001, x$degreesoffreedom)
   stop  <- qt(0.999, x$degreesoffreedom)
   xx    <- seq(start, stop, length=1001)
   if (x$P==0) {
      plot(xx, rep(1,1001), xlab = "", ylab="", type="l")
   } else {
      y     <- dt(xx, x$degreesoffreedom)
      xx    <- xx/sqrt(x$P) + x$expectation
      plot(xx, y, xlab = "", ylab="", type="l")
   }
}

summary.tdistribution <- function(object, ...) {
print(object)
if (object$degreesoffreedom>2) 
   cat("Variance ", object$degreesoffreedom/(object$degreesoffreedom-2)/object$P, "\n")
else
   cat("This distribution has no variance.\n")
cat("Mode ", object$expectation, "\n")
}

simulate.tdistribution <- function(object, nsim = 1, ...) {
  res <- rt(nsim, object$degreesoffreedom)
  res/sqrt(object$P) + object$expectation
}

difference.tdistribution <- function(object1, object2) {
if (class(object2)[1]=="normal") {
    cat("WARNING: The resulting distribution is an approximation.\n")
    tdistribution(object2$expectation - object1$expectation, 
    object1$degreesoffreedom*(1+object1$P/object2$precision)^2, 
    P = 1/(1/object1$P + 1/object2/precision))
} else if (class(object2)[1]=="tdistribution") {
    cat("WARNING: The resulting distribution is an approximation.\n")
    tdistribution(object2$expectation - object1$expectation,
    (1/object1$P + 1/object2$P)^2/
    (1/object1$P^2/object1$degreesoffreedom + 1/object2$P^2/object2$degreesoffreedom),
    P = 1/(1/object1$P + 1/object2$P))
} else
cat("ERROR: The difference function is not implemented for this pair of objects.\n")
}

credibilityinterval.tdistribution <- function(object, prob=0.95) {
  res <- qt((1-prob)/2, object$degreesoffreedom)
  c(res, -res)/sqrt(object$P) + object$expectation
}

p.value.tdistribution <- function(object, point=0) {
# For any distribution, this should return the probability outside the 
# credibility region with boundary containing the origin. 
  pt((point - object$expectation)*sqrt(object$P), object$degreesoffreedom, 
  	  lower.tail = (point < object$expectation))*2
}


##################################
## CLASS mtdistribution 
##################################
# Data: 
# expectation: must be of length k>=1
# OPTIONAL: names: character vector of length k, giving the names of the dimensions
# degreesoffreedom
# P 

mtdistribution <- function(expectation=c(0,0), degreesoffreedom = 10000, P = diag(length(expectation))) {
if (degreesoffreedom < 1) 
  cat("ERROR: The degrees of freedom must be at least 1.\n") 
else if (any(eigen(P)$values<0))
  cat("ERROR: The P matrix must be symmetric and non-negative definite.\n")
else {
result <- list(expectation=expectation, degreesoffreedom=degreesoffreedom, P=P)
class(result) <- c("mtdistribution", "probabilitydistribution")
result
}
}

expectation.mtdistribution <- function(object) {
   object$expectation
}

variance.mtdistribution <- function(object) {
   result <- rep(object$degreesoffreedom/(object$degreesoffreedom-2), prod(dim(object$P)))*ginv(object$P)
   0.5*(result + t(result))
}   

precision.mtdistribution <- function(object) {
   rep((object$degreesoffreedom-2)/object$degreesoffreedom, prod(dim(object$P)))*object$P
}

marginal.mtdistribution <- function(object, v) {
# v gives the indices for the dimensions which should be kept
   k <- length(object$expectation)
   v <- as.integer(v)
   if(any(duplicated(v)))
      cat("ERROR: Duplicated indices not allowed.\n")
   else if (max(v)>k)
      cat("ERROR: Too large index.\n")
   else if (min(v)<1)
      cat("ERROR: Too small index.\n")
   else {
      object$expectation <- object$expectation[v]
      if (length(v)==k)
         object$P  <- object$P[v,v]
      else {
         object$P  <- object$P[v,v] - object$P[v,-v]%*%ginv(object$P[-v,-v])%*%object$P[-v,v]  
         object$P <- 0.5*(object$P + t(object$P))
      }    
      if (length(object$expectation)==1)
         tdistribution(object$expectation, object$degreesoffreedom, P = object$P)
      else 
         object
   }
}

conditional.mtdistribution <- function(object, v, val) {
# v gives the indices for the dimensions whose values should be fixed, and the 
# fixed values are given by val.
   if (missing(v)) 
      cat("ERROR: Argument giving indices for conditioning is missing.\n")
   else if (missing(val))
      cat("ERROR: Argument giving values for conditioning is missing.\n")
   else { 
   k <- length(object$expectation)
   v <- as.integer(v)
   if (length(v)<1) 
      cat("ERROR: Cannot condition on empty index vector.\n")
   else if (any(duplicated(v)))
      cat("ERROR: Duplicated indices not allowed.\n")
   else if (max(v)>k)
      cat("ERROR: Too large index.\n")
   else if (min(v)<1)
      cat("ERROR: Too small index.\n")
   else if (length(v)==k)
      cat("ERROR: Cannot condition on all indices simultaneously.\n")
   else if (length(val)!=length(v))
      cat("ERROR: The number of values must be equal to the number of indices.\n")
   else {
      newexp <- object$expectation[-v] - ginv(object$P[-v,-v])%*%object$P[-v,v]%*%(val - object[v])
      newPfactor <- object$degreesoffreedom/(object$degreesoffreedom + 
         (val - object$expectation[v])%*%
         (object$P[v,v]-object$P[v,-v]%*%ginv(object$P[-v,-v])%*%object$P[-v,v])%*%
         (val - object$expectation[v]))
      object$P <- object$P[-v,-v]
      object$P <- rep(newPfactor, prod(dim(object$P)))*object$P
      object$P <- 0.5*(object$P + t(object$P))
      object$expectation <- newexp
      if (length(object$expectation)==1)
         class(object) <- "tdistribution"
      object
   }
   }
}

probabilitydensity.mtdistribution <- function(object, val, log=FALSE, normalize=TRUE) {
    k <- length(object$expectation)
    val <- as.vector(val)
    logintconstant <- lgamma(object$degreesoffreedom/2) - lgamma((object$degreesoffreedom + 
        k)/2) +  k/2*log(pi) -0.5*log(det(object$P)) - object$degreesoffreedom/2*log(object$degreesoffreedom)
    logdensity <- -(object$degreesoffreedom + k)/2*log(object$degreesoffreedom + 
        (val - object$expectation)%*%object$P%*%(val-object$expectation))
    ifelse(log, logdensity - logintconstant, exp(logdensity - logintconstant))
}

print.mtdistribution <- function(x, ...) {
cat("A Multivariate t probability distribution.\n")
cat("Expectation ", x$expectation, "\n")
cat("Degrees of freedom ", x$degreesoffreedom, "\n")
cat("Parameter proportional to precision ", x$P, "\n")
}

plot.mtdistribution <- function(x, onlybivariate=FALSE, ...) {
   if (onlybivariate & length(x$expectation)==2) {
      # To plot improper distributions reasonably: 
      Ptmp <- eigen(x$P)
      v <- ifelse(Ptmp$values, Ptmp$values, 1)
      Pnew <- Ptmp$vectors%*%diag(v)%*%t(Ptmp$vectors)
      st1 <- sqrt(Pnew[1,1]-Pnew[1,2]*Pnew[2,1]/Pnew[2,2])
      st2 <- sqrt(Pnew[2,2]-Pnew[2,1]*Pnew[1,2]/Pnew[1,1])
      fxt    <- qt(0.999, x$degreesoffreedom)
      startA <- x$expectation[1] - fxt/st1
      stopA  <- x$expectation[1] + fxt/st1
      startB <- x$expectation[2] - fxt/st2
      stopB  <- x$expectation[2] + fxt/st2
#
#      A <- marginal(x, 1)
#      B <- marginal(x, 2)
#      startA <- qt(0.001, A$degreesoffreedom)/sqrt(A$P)+A$expectation
#      stopA  <- qt(0.999, A$degreesoffreedom)/sqrt(A$P)+A$expectation
#      startB <- qt(0.001, B$degreesoffreedom)/sqrt(B$P)+B$expectation
#      stopB  <- qt(0.999, B$degreesoffreedom)/sqrt(B$P)+B$expectation
      logf <- function(xx, ...) {
         -(x$degreesoffreedom + 2)/2*log(x$degreesoffreedom + (xx-x$expectation) %*%
                 x$P %*% (xx-x$expectation))
      }
      densitycontour(logf, c(startA, stopA, startB, stopB), data=0)
   } else {
      k <- length(x$expectation)
      old.par <- par(no.readonly = TRUE); on.exit(par(old.par))
      par(mfcol=c(k,k))
      for (i in 1:k) for (j in 1:k) {
         if (i==j) 
            plot(marginal(x, i))
         else 
            plot(marginal(x, c(i,j)), TRUE)
      }
   }
}

summary.mtdistribution <- function(object, ...) {
print(object)
}

simulate.mtdistribution <- function(object, nsim = 1, ...) {
   helper <- rgamma(nsim, object$degreesoffreedom/2, object$degreesoffreedom/2)
   k <- length(object$expectation)
   res <- simulate(mnormal(rep(0, k), P=object$P), nsim)
   res*matrix(rep(1/sqrt(helper), k), nsim, k) + matrix(object$expectation, nsim, k, byrow=TRUE)
}

p.value.mtdistribution <- function(object, point=0) {
# For any distribution, this should return the probability outside the 
# credibility region with boundary containing the origin. 
    k <- length(object$expectation)
    if (sum(point*point)==0) point <- rep(0, length(object$expectation))
    1-pf(as.vector((object$expectation-point)%*%object$P%*%(object$expectation-point))/k, k, object$degreesoffreedom)
}

##################################
## CLASS fdistribution 
##################################
# Data: 
# Parameters: df1 and df2

fdistribution <- function(df1 = 1, df2 = 1) {
if (df1<0) 
cat("ERROR in specification\n")
else if (df2<0)
cat("ERROR in specification\n")
else {
	result <- list(df1=df1, df2=df2)
	class(result) <- c("fdistribution", "probabilitydistribution")
	result
}
}

expectation.fdistribution <- function(object) {
if (object$df2<=2) 
   cat("ERROR: The expectation does not exist.\n")
else
   object$df2/(object$df2-2)
}

variance.fdistribution <- function(object) {
if (object$df2<=4)
   cat("ERROR: The variance does not exist.\n")
else 
   2*object$df2^2*(object$df1+object$df2-2)/
   (object$df1*(object$df2-2)^2*(object$df2-4))
}

precision.fdistribution <- function(object) {
if (object$df2<=4)
   cat("ERROR: The variance does not exist.\n")
else 
   1/(2*object$df2^2*(object$df1+object$df2-2)/
   (object$df1*(object$df2-2)^2*(object$df2-4)))
}

probabilitydensity.fdistribution <- function(object, val, log=FALSE, normalize=TRUE) {
   df(val, object$df1, object$df2, log=log)
}

cdf.fdistribution <- function(object, val) {
   pf(val, object$df1, object$df2)
}

invcdf.fdistribution <- function(object, val) {
   if (val<=0 | val>=1) 
      cat("ERROR: Illegal input for the invcdf function.\n")
   else 
      qf(val, object$df1, object$df2)
}

print.fdistribution <- function(x, ...) {
cat("An F probability distribution.\n")
cat("First degre of freedom: ", x$df1, "\n")
cat("Second degre of freedom: ", x$df2, "\n")
}

plot.fdistribution <- function(x, ...) {
   start <- qf(0.0001, x$df1, x$df2)
   stop  <- qf(0.999, x$df1, x$df2)
   xx    <- seq(start, stop, length=1001)
   yy    <- df(xx, x$df1, x$df2)
   plot(xx, yy, xlab = "", ylab="", type="l")
}

summary.fdistribution <- function(object, ...) {
print(object)
cat("Expectation: ", expectation(object), "\n")
cat("Variance: ", variance(object), "\n")
}

simulate.fdistribution <- function(object, nsim=1, ...) {
   rf(nsim, object$df1, object$df2)
}

p.value.fdistribution <- function(object, point=0) {
# For any distribution, this should return the probability outside the 
# credibility region with boundary containing the point. 
    res <- pf(point, object$df1, object$df2) 
    if (res > 0.5) res <- 1-res
    2*res
}

credibilityinterval.fdistribution <- function(object, prob=0.95) {
c(qf((1-prob)/2, object$df1, object$df2), 
  qf(1-(1-prob)/2, object$df1, object$df2))
}


##################################
## CLASS wishart
##################################
# Data: 
# S: symmetric non-negative definite matrix. 
# df: degrees of freedom

wishart <- function(S = diag(2), df = 1) {
   if (any(t(S)!=S)) 
      cat("ERROR: The S matrix must be symmetric.\n")
   else if (min(eigen(S)$values)<0)
      cat("ERROR: The S matrix must be non-negative definite.\n")
   else if (df < 0) 
      cat("ERROR: The degrees of freedom cannot be negative.\n")
   else {
      result <- list(S = S, df = df)
      class(result) <- c("wishart", "probabilitydistribution")
      result
   }
}

expectation.wishart <- function(object) {
   object$df*ginv(object$S)
}

variance.wishart <- function(object) {
   cat("Not implemented yet.\n")
}

precision.wishart <- function(object) {
   cat("Not implemented yet.\n")
}

marginal.wishart <- function(object, v) { 
   # v gives the indices for the dimensions which should be kept. 
   k <- dim(object$S)[1]
   v <- as.integer(v)
   if (any(duplicated(v)))
      cat("ERROR: Duplicated indices not allowed.\n")
   else if (max(v)>k)
      cat("ERROR: Too large index.\n")
   else if (min(v)<1)
      cat("ERROR: Too small index.\n")
   else if (length(v)==k) {
      # permute the indices
      object$S <- object$S[v,v]
      object
   }
   else {
      object$S   <- object$S[v,v] - object$S[v,-v] %*% ginv(object$S[-v,-v]) %*% object$S[-v, v]
          object$S <- 0.5*(object$S + t(object$S))
      if (length(object$S)==1)
         gammadistribution(object$df/2, object$S/2)
      else 
         object
   }
}

conditional.wishart <- function(object, v, val) {
# v gives the indices for the dimensions whose values should be fixed, and the 
# fixed values are given by val. 
   cat("Not implemented yet.\n")
}

probabilitydensity.wishart <- function(object, val, log=FALSE, normalize=TRUE) {
   k <- dim(object$S)[1]
   result <- det(val)^((object$df - k - 1)/2)*exp(-sum(diag(object$S%*%val))/2)
   if (normalize) 
      result <- result/(2^(object$df*k/2)*pi^(k*(k-1)/4)*prod(gamma((object*df + 1 - (1:k))/2)))
   ifelse(log, log(result), result)
}

print.wishart <- function(x, ...) {
cat("A Wishart probability distribution.\n")
cat("Degrees of freedom: ", x$df, "\n")
cat("S matrix:\n")
print(x$S) 
}

plot.wishart <- function(x, onlybivariate=FALSE, ...) {
   if (onlybivariate & length(x$expectation)==2) {
      # To plot improper distributions reasonably: 
      # we need to fill in code here. For now, assume that 
      # S is positive definite. 
      k <- dim(x$S)[1]
      A <- marginal(x, 1)
      B <- marginal(x, 2)
      startA <- invcdf(A, 0.001)
      stopA  <- invcdf(A, 0.999)
      startB <- invcdf(B, 0.001)
      stopB  <- invcdf(B, 0.999)
      logf <- function(X, ...) {
         0.5*(x$df - k - 1)*log(det(X)) - 0.5*sum(diag(x$S%*%X))
      }
      densitycontour(logf, c(startA, stopA, startB, stopB), data=0)
   } else {
      k <- dim(x$S)[1]
      old.par <- par(no.readonly = TRUE); on.exit(par(old.par))
      par(mfcol=c(k,k))
      for (i in 1:k) for (j in 1:k) 
         if (i==j) plot(marginal(x, i))
         else plot(marginal(x, c(i,j)), TRUE)
   }
}

summary.wishart <- function(object, ...) {
print(object)
}




#To do: BELOW


simulate.wishart <- function(object, nsim = 1, ...) {
	if (det(object$P)==0) {
		cat("ERROR: Cannot simulate from improper distribution")
	} else {
		k <- length(object$expectation)
		B <- matrix(rnorm(nsim*k), nsim, k)%*%chol(solve(object$P))
		B + matrix(object$expectation, nsim, k, byrow=TRUE)
	}
}

linearpredict.wishart <- function(object, X = rep(1,length(object$expectation)), PP = diag(length(X)/length(object$expectation)), ...) {
   if (is.matrix(X) && dim(X)[2] != length(object$expectation))
      cat("ERROR: X has wrong number of columns.\n")
   else if (!is.matrix(X) && length(X) != length(object$expectation))
      cat("ERROR: X has wrong length.\n")
   else if (dim(PP)[1] != length(X)/length(object$expectation) | dim(PP)[2] != length(X)/length(object$expectation))
      cat("ERROR: P as wrong format.\n")
   else {
      k <- length(object$expectation)
      s <- length(X)/k
      X <- matrix(X, s, k)
      expectation <- c(X%*%object$expectation, object$expectation)
      P   <- matrix(NA, s+k, s+k)
      P[1:s,1:s] <- PP
      P[1:s,s+(1:k)] <- -PP%*%X
      P[s+(1:k),1:s] <- -t(PP%*%X)
      P[s+(1:k),s+(1:k)] <- t(X)%*%PP%*%X + object$P
      result <- list(expectation=expectation, P=P)
      class(result) <- c("mnormal", "probabilitydistribution")
      result
   }
}









####################################################################
####################################################################
## EXTRA FUNCTIONS: 

designOneGroup <- function(n) {
matrix(rep(1,n),n,1)
}

designTwoGroups <- function(n,m) {
matrix(c(rep(1,n+m), rep(0,n), rep(1,m)), n+m, 2)
}

designManyGroups <- function(v) {
result <- matrix(0, sum(v), length(v))
w <- cumsum(v)
result[,1] <- 1
for (i in 2:length(v)) result[(w[i-1]+1):w[i], i] <- 1
result
}

designBalanced <- function(factors, replications = 1, interactions = FALSE) {
# factors is a vector with length equal to the number of dimensions (or factors), 
# and where each value is the number of levels for each factor. 
# replications is the number of replications for each combination of factors
# interactions can have many forms: If it has length 1 it is either true or false
# MORE FORMS WILL BE IMPLEMENTED LATER?? but users can always delete columns they don't like.
# Make the design without interaction first, then, possibly add the interactions
if (min(factors)<2) 
   cat("All factor level counts must be at least 2.\n")
else {
   m <- prod(factors)*replications
   n <- length(factors) 
   result <- matrix(0, m, sum(factors - 1) + 1)
   result[,1] <- 1
   counter <- 2
   rep1 <- m
   for (i in 1:n) {
      rep1 <- rep1/factors[i]
      for (j in 2:factors[i]) {
         r1 <- rep(0, rep1*factors[i])
	 r1[rep1*(j-1)+(1:rep1)] <- 1
         result[,counter] <- r1
         counter <- counter + 1
      }
   }
   if (interactions) {
      csum <- c(0, cumsum(factors-1))
      result <- cbind(result, matrix(0, m, prod(factors) - sum(factors - 1) - 1))
      ccount <- rep(1, n)
      repeat {
         if (sum(ccount>1)>1) {
            result[,counter] <- result[,1]
            for (i in 1:n) if (ccount[i]>1) result[,counter] <- result[,counter]*result[,csum[i]+ccount[i]]
	    counter <- counter + 1
         }
	 if (sum(ccount)==sum(factors)) break
         k <- 1
         while (ccount[k] == factors[k]) k <- k + 1
         ccount[k] <- ccount[k]+1
	 if (k>1) ccount[1:(k-1)] <- 1
      } 
   }
   result 
   }
}


designFactorial <- function(nfactors, replications = 1, interactions = FALSE) {
# Creates a design matrix for two-level factors. 
   n <- 2^nfactors*replications
   result <- matrix(1, n, nfactors + 1)
   for (i in 1:nfactors) {
      res <- c(rep(-1, 2^(nfactors-i)*replications), rep(1, 2^(nfactors-i)*replications))
      result[,i+1] <- res
   }
   if (interactions) {
      counter <- nfactors + 2
      result <- cbind(result, matrix(1, n, 2^nfactors - nfactors - 1))
      for (cc in 2:nfactors) {
         #go through all ways to select cc factors among nfactors to participate in the interaction
         selection <- 1:cc
         repeat {
            result[,counter] <- apply(result[,selection + 1], 1, prod)
            counter <- counter + 1
            if (selection[1] == nfactors - cc + 1) break
            j <- cc
            while (selection[j] == nfactors - cc + j) j <- j - 1
            selection[j] <- selection[j] + 1
            if (j < cc) selection[(j+1):cc] <- selection[j] + (1:(cc-j))
         }
      }
   }
   result 
}

### NICE TO HAVE: 
#A function that takes a set of factors (R standard definition) 
#and produces (suggests) a design matrix (with column names) (doesn't R 
#do this already?) 

# Should be contained in the package description: A discussion of how these 
# functions relate to the standard lm functions in R. 

######################################################################
# Special functions: 

flat <- function(k) {
# Produces a flat uniform distribution of dimension k
if (k==1) 
   uniformdistribution(-Inf, Inf)
else 
   muniformdistribution(rep(-Inf, k), rep(Inf, k))
}

linearmodel <- function(data, design) {
if (length(data)!=dim(design)[1]) 
   cat("ERROR: The design matrix has wrong format.\n")
else {
result <- list()
result$P <- t(design)%*%design
result$mu <- ginv(result$P) %*% t(design) %*% data
n <- dim(design)[1]
k <- dim(design)[2]
result$alpha <- (n-k)/2
result$beta  <- (sum(data^2) -  data %*% design %*% result$mu)/2
if (k>1) 
   class(result) <- c("mnormalexpgamma", "probabilitydistribution") 
else
   class(result) <- c("normalexpgamma", "probabilitydistribution")
result
}
}

posteriornormal1 <- function(data) {
   linearmodel(data, designOneGroup(length(data)))
}

posteriornormal2 <- function(data1, data2) {
   marginal(linearmodel(c(data1, data2), designTwoGroups(length(data1), length(data2))), 2:3)
}

anovatable <- function(data, design, subdivisions=c(1, dim(design)[2]-1)) {
if (sum(subdivisions)!=dim(design)[2]) 
cat("ERROR: The subdivision vector is not correct.\n")
else if (subdivisions[1] != 1) 
cat("ERROR: The subdivision vector is not correct.\n")
else {
k <- length(subdivisions)
kk <- cumsum(subdivisions)
result <- matrix(NA, k+1, 5)
for (i in 2:k) {
   result[i-1,1] <- sum((fittedvalues(data, design[,1:kk[i-1], drop=FALSE]) - 
     fittedvalues(data, design[, 1:kk[i], drop=FALSE]))^2)
   result[i-1,2] <- subdivisions[i]
}
result[k, 1] <- sum((data - fittedvalues(data, design))^2)
result[k, 2] <- dim(design)[1]-dim(design)[2]
result[k+1, 1] <- sum(result[1:k,1])
result[k+1, 2] <- dim(design)[1] - 1
result[1:k,3] <- result[1:k,1]/result[1:k,2]
result[1:(k-1),4] <- result[1:(k-1),3]/result[k,3]
result[1:(k-1),5] <- pf(result[1:(k-1),4], result[1:(k-1),2], result[k,2], lower.tail=FALSE)
result
}
}


leastsquares <- function(data, design) {
if (length(data)!=dim(design)[1]) 
   cat("ERROR: The design matrix has wrong format.\n")
else 
   as.vector(ginv(t(design)%*%design)%*%t(design)%*%data)
}

fittedvalues <- function(data, design) {
if (length(data)!=dim(design)[1]) 
   cat("ERROR: The design matrix has wrong format.\n")
else 
   as.vector(design%*%ginv(t(design)%*%design)%*%t(design)%*%data)
}

densitycontour <- function (logf, limits, data, ...) 
{
# This function is adapted from mycontour in the LearnBayes package. 
    LOGF = function(theta, data) {
        if (is.matrix(theta) == TRUE) {
            val = matrix(0, c(dim(theta)[1], 1))
            for (j in 1:dim(theta)[1]) val[j] = logf(theta[j, 
                ], data)
        }
        else val = logf(theta, data)
        return(val)
    }
    ng = 50
    x0 = seq(limits[1], limits[2], len = ng)
    y0 = seq(limits[3], limits[4], len = ng)
    X = outer(x0, rep(1, ng))
    Y = outer(rep(1, ng), y0)
    n2 = ng^2
    Z = LOGF(cbind(X[1:n2], Y[1:n2]), data)
    Z = Z - max(Z)
    Z = matrix(Z, c(ng, ng))
    contour(x0, y0, Z, levels = seq(-6.9, 0, by = 2.3), drawlabels=FALSE, lwd = 2, 
        ...)
}


######################################################################

# Functions that are REMOVED, for now: 
#lineartransform <- function(object, ...) {
#UseMethod("lineartransform")
#}
#lineartransform.default <- function(object, ...) {
#print("The lineartransform function is not implemented for this object.")
#}
#lineartransform.normal <- function(object, v) {
#if (length(v)>1) print("ERROR: Contrast vector of incorrect length.")
#else if (v==0) print("ERROR: Cannot deal with objects with infinite precision.") 
#else {
#object$expectation <- object$expectation * v
#object$precision   <- object$precision / v^2
#object  
#}
#}
#lineartransform.mnormal <- function(object, A) {
#	k <- length(object$expectation)
#	if (k==length(A)) {
#		object$expectation <- sum(object$expectation*A)
#		object$precision <- 1/(A%*%solve(object$precision)%*%A)
#		class(object)[1] <- "normal"
#		object
#	} else if (k==dim(A)[1]) {
#		object$expectation <- object$expectation%*%A
#		object$precision <- solve(t(A)%*%solve(object$precision)%*%A)
#		object		
#	} else {
#		 cat("ERROR: Transformation of incorrect dimension.\n")
#	}
#}
#lineartransform.gammadistribution <- function(object, v) {
#if (length(v)>1) 
#   cat("ERROR: Contrast vector of incorrect length.\n")
#else if (v==0) 
#   cat("ERROR: Cannot deal with objects with infinite precision.\n") 
#else {
#   object$beta <- object$beta/v
#   object
#}
#}
#lineartransform.tdistribution <- function(object, v) {
#if (length(v)>1) print("ERROR: Contrast vector of incorrect length.")
#else if (v==0) print("ERROR: Cannot deal with objects with infinite precision.") 
#else {
#object$expectation <- object$expectation * v
#object$P <- object$P/sqrt(v)
#object  
#}
#}

#quotient <- function(object1, object2,  ...) {
# This should mainly be used to generate the F-distribution?? that appears when 
# studying the quotient of two independent chi-square distributions. 
#UseMethod("quotient")
#}
#quotient.default <- function(object, ...) {
#print("The quotient function is not implemented for this object.")
#}

