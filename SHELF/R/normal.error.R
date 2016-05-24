normal.error <-
function(parameters, values, probabilities, weights){
	sum(weights * (pnorm(values, parameters[1], exp(parameters[2])) - probabilities)^2)
}
