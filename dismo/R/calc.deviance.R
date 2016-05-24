# j. leathwick/j. elith
#
# version 2.1 - 5th Sept 2005
#
# function to calculate deviance given two vectors of raw and fitted values
# requires a family argument which is set to binomial by default
#
#

calc.deviance <-  function(obs, pred, weights = rep(1,length(obs)), family="binomial", calc.mean = TRUE) {

if (length(obs) != length(pred)) {   stop("observations and predictions must be of equal length") }

y_i <- obs
u_i <- pred
 
family = tolower(family)
 
if (family == "binomial" | family == "bernoulli") {
 
   deviance.contribs <- (y_i * log(u_i)) + ((1-y_i) * log(1 - u_i))
   deviance <- -2 * sum(deviance.contribs * weights)

} else if (family == "poisson") {

    deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - (y_i - u_i)
    deviance <- 2 * sum(deviance.contribs * weights)

} else if (family == "laplace") {

    deviance <- sum(abs(y_i - u_i))
	
} else if (family == "gaussian") {

    deviance <- sum((y_i - u_i) * (y_i - u_i))
	
} else {
	stop('unknown family, should be one of: "binomial", "bernoulli", "poisson", "laplace", "gaussian"')
}

if (calc.mean) deviance <- deviance/length(obs)

return(deviance)

}

