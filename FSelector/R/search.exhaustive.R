exhaustive.search <- function(attributes, eval.fun) {
	len = length(attributes)
	if(len == 0)
		stop("Attributes not specified")
	
	eval.fun = match.fun(eval.fun)
	best = list(
		result = -Inf,
		attrs = rep(0, len)
	)
	
	# main loop
	# for each subset size
	for(size in 1:len) {
		child_comb = combn(1:len, size)
		# for each child
		for(i in 1:dim(child_comb)[2]) {
			subset = rep(0, len)
			subset[child_comb[, i]] = 1
			result = eval.fun(attributes[as.logical(subset)])
			if(result > best$result) {
				best$result = result
				best$attrs = subset
			}
		}
	}
	return(attributes[as.logical(best$attrs)])
}