make.thresholds.matrix <- function(item.params, design.matrix = "normal", make.from = "deltas", theta.interval = c(-10, 10), throld = 0.5, 
	alpha = 1, c.params = 0, ...) {
	#print("default")
	
	# Provides a predicted probability of a given response for a polytomous
	# item. Mostly used in other functions.

	predicted.prob = function(theta, response, design.matrix, parameters, slope) {
		numerator =  exp(slope * (theta * (response - 1) + sum(parameters * design.matrix[response, ])))
		denominator = 0
		for (k in 1:nrow(design.matrix)) {
			denominator = denominator + exp(slope * (theta * (k - 1) + sum(parameters * design.matrix[k, ])))
		}
		return(numerator/denominator)
	}

	throlds.minimize = function(theta, threshold, slope) {
		numerator = exp(slope * (theta - threshold))
		denominator = 1 + exp(slope * (theta - threshold))
		prob <- numerator/denominator
		return((prob - throld)^2)
	}

	# The function passed to R's 'optimize' function. Finds the cumulative
	# predicted probability of a response, then returns (predict P - .5)^2.
# Minimizing this function finds the Thurstone threshold.
minimize.fun = function(theta, response, design.matrix, parameters, slope, guess) {
		max.score = nrow(design.matrix)
		categories = max.score:response
		predictions = rep(NA, times = max.score)
		for (x in categories) {
			predictions[x] = guess + (1 - guess) * predicted.prob(theta, x, design.matrix, parameters, slope)
		}
		total = sum(predictions, na.rm = TRUE)
		sq.err = (total - throld)^2
		return(sq.err)
	}

	# First creates an item design matrix if one not provided. Options
	# are "normal", "conquest", and an actual matrix. Then minimizes
# 'minimize.fun' for each threshold for the item.
get.thresholds = function(parameters, design.matrix = "normal", theta.interval = c(-10, 10), slope, guess) {
		#print(parameters)
		#print("get")
		parameters = as.numeric(parameters)
		max.length <- length(parameters)
		parameters = parameters[!is.na(parameters)]
		n.parameters = length(parameters)
		thresholds = rep(NA, times = max.length)
		if(make.from == "deltas") {
		
		n.categories = n.parameters + 1

		if (is.character(design.matrix)) {
			if (design.matrix == "conquest") {
				design.matrix = matrix(0, ncol = n.parameters, nrow = n.categories)
				design.matrix[, 1] = -1 * (1:n.categories - 1)
				i = 2
				while (i <= n.categories - 1) {
					design.matrix[i, 2:i] = -1
					i = i + 1
				}
			} else {
				design.matrix = matrix(0, ncol = n.parameters, nrow = n.categories)
				i = 2
				while (i <= n.categories) {
					design.matrix[i, 1:(i - 1)] = -1
					i = i + 1
				}
			}
			#print(design.matrix)
			}


		for (response in 2:n.categories) {
			opti = optimize(minimize.fun, interval = theta.interval, response = response, design.matrix = design.matrix, parameters = parameters, 
				slope = slope, guess = guess)
			thresholds[response - 1] = opti$minimum
			if (opti$objective > 10^(-5)) 
				message(paste("objective = ", opti$objective))
		}
		}
		else if(make.from == "thresholds"){
			for (response in 1:n.parameters) {
				opti <- optimize(throlds.minimize,interval = theta.interval, threshold = parameters[response], slope = slope)
				thresholds[response] <- opti$minimum
				if (opti$objective > 10^(-5)) 
				message(paste("objective = ", opti$objective))
				
			}
		}

		return(thresholds)

	}

	# Accepts a matrix of item parameters (rows are items) and feeds the items
	# to 'get.thresholds'. NAs are fine in the parameter.matrix, allowing
# items to have differing numbers of categories. design.matrix may be a
# single choice for all items or a list of matrices, one for each item.
apply.thresholds = function(parameter.matrix, design.matrix = "normal", theta.interval = c(-10, 10), slope = slope, guess = guess) {


		threshold.matrix <- as.matrix(mapply(get.thresholds, as.data.frame(t(parameter.matrix)), t(design.matrix), slope = slope, guess = guess))
		if (NCOL(threshold.matrix) > 1) 
			threshold.matrix <- t(threshold.matrix)

		return(threshold.matrix)

	}
	if(any(c.params != 0)) {
		if(make.from != "deltas")
			stop("Cannot run 3PL using threshold inputs")
		if(NCOL(item.params) != 1)
			stop("Cannot run polytomous 3PL")
	}
	if(make.from == "deltas")
		message("Assuming partial credit model")
	else if(make.from == "thresholds")
		message("Assuming graded response model")

	return(apply.thresholds(item.params, design.matrix, theta.interval, alpha, c.params))


}
