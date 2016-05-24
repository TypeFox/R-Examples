
forward.search <- function(attributes, eval.fun) {
	return(greedy.search(attributes, eval.fun, TRUE))
}

backward.search <- function(attributes, eval.fun) {
	return(greedy.search(attributes, eval.fun, FALSE))
}

greedy.search <- function(attributes, eval.fun, forward = TRUE) {
	if(length(attributes) == 0)
		stop("Attributes not specified")
	
	eval.fun = match.fun(eval.fun)
	best = list(
		result = -Inf,
		attrs = rep(as.numeric(!forward), length(attributes))
	)
	
	# initial evaluation for full set when backward
	if(!forward) {
		best$result = eval.fun(attributes[as.logical(best$attrs)])
	}
	
	forward_text = ifelse(forward, "forward", "backward")
	
	# main loop
	repeat {
		# create new matrix of children to evaluate
		children = create.children(best$attrs, forward_text)
		if(is.null(children))
			break()
		
		# evaluate and find the best of them
		children_results = apply(children, 1, function(vec) {
				eval.fun(attributes[as.logical(vec)])
			})
		local_best = find.best(children_results)
		
		# compare to the best so far
		if(local_best$result > best$result) {
			best$result = local_best$result
			best$attrs = children[local_best$idx,]
		} else {
			break()
		}
	}
	
	return(attributes[as.logical(best$attrs)])

}
