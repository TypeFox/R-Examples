t.error <-
function(parameters, values, probabilities, weights, degreesfreedom){
	sum(weights * (pt((values-parameters[1]) / exp(parameters[2]), degreesfreedom) - probabilities)^2)
}
