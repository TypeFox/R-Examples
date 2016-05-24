gamma.error <-
function(parameters, values, probabilities, weights){
	sum(weights * (pgamma(values, exp(parameters[1]), exp(parameters[2])) -probabilities)^2)
}
