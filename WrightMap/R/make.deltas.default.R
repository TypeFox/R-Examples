make.deltas.default <-
function(item.params, cross.params = 0, step.params = 0, item.sign = 1, step.sign = 1, inter.sign = 1, ...) {
	#print(item.params)
	#print(step.params)
throlds <- item.params
	crosses <- cross.params
	steps <- step.params
	#print(throlds)
	#print(crosses)
	num.items <- length(throlds)
	if (length(steps) > 1) 
		steps <- matrix(steps, nrow = num.items, ncol = length(steps), byrow = TRUE)
	#print(steps)
	
	if (length(crosses) > 1) {
		item.nums <- 1:num.items
		full.steps <- matrix(nrow = max(item.nums), ncol = ncol(crosses))
		full.steps[, 2] <- 0

		full.steps[item.nums %in% unlist(crosses[1]), ] <- as.matrix(crosses)


		crosses <- full.steps[, -1]
	} else crosses <- 0
	#print(inter.sign)
	#print(item.sign)
	#print(step.sign)
	throlds <- crosses * inter.sign + throlds * item.sign + steps * step.sign
	throlds <- throlds[rowSums(!is.na(throlds)) != 0, ]

	return(throlds)


}
