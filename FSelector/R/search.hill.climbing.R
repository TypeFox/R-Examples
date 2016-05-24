hill.climbing.search <- function(attributes, eval.fun) {
	if(length(attributes) == 0)
		stop("Attributes not specified")
	
	eval.fun = match.fun(eval.fun)
	best = list(
		result = -Inf,
		attrs = rep(0, length(attributes))
	)
	while(sum(best$attrs) == 0)
		best$attrs = sample(c(0,1), length(attributes), replace = TRUE)
	best$result = eval.fun(attributes[as.logical(best$attrs)])
	
	evaluated_states = list(
		attrs = matrix(best$attrs, nrow = 1, ncol = length(attributes), byrow = TRUE),
		results = best$result
	)
	
	eval_state <- function(state, evaluated_states) {
		idx = find.subset(evaluated_states$attrs, state)
		if(length(idx) == 0) {  # needs to be evaluated
			return(list(
				to_be_saved = TRUE,
				result = eval.fun(attributes[as.logical(state)])
				))
		} else if(length(idx) == 1) { # already evaluated
			return(list(
				to_be_saved = FALSE,
				result = evaluated_states$results[idx]
				))
		} else {
			stop("Internal error")
		}
	}
	
	# main loop
	repeat {
		# find neighbours
		children = create.children(best$attrs, direction = "both")
		if(is.null(children))
			break()
		
		# evaluate and find the best of them
		children_evaluated = apply(children, 1, function(vec) {
				eval_state(vec, evaluated_states)
			})

		children_results = sapply(children_evaluated, function(x) x$result)
		children_to_be_saved = sapply(children_evaluated, function(x) x$to_be_saved)
		
		# save children at evaluated_states
		evaluated_states$attrs = rbind(evaluated_states$attrs, children[children_to_be_saved,, drop = FALSE])
		evaluated_states$results = c(evaluated_states$results, children_results[children_to_be_saved])
		
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