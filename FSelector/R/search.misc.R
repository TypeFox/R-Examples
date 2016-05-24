# returns indicies
find.subset <- function(subsets.matrix, subset) {
	subset = as.vector(subset)
	len = length(subset)
	if(len == 0)
		stop("Empty atrributes subset.")
	if(dim(subsets.matrix)[2] != len)
		stop("Incorrect dimensions.")
		
	if(dim(subsets.matrix)[1] == 0)
		return(as.integer(NULL))
	
	cond = rep(TRUE, dim(subsets.matrix)[1])
	for(i in 1:len)
		cond = cond & (subsets.matrix[,i] == subset[i])
	return(which(cond))
}


create.children <- function(parent, direction = c("forward", "backward", "both"), omit.func = NULL ) {
	direction = match.arg(direction)

	if(!is.null(omit.func)) {
		omit.func = match.fun(omit.func)
	}

	cols = length(parent)
	if(cols <= 0)
		stop("Parent attribute set cannot be empty.")
	
	m1 = NULL
	m2 = NULL
	
	if(direction == "forward" || direction == "both") {
		rows = cols - sum(parent)
		if(rows > 0) {
			m1 = matrix(parent, ncol = cols, nrow = rows, byrow = TRUE)

			current_row = 1
			current_col = 1
			repeat {
				if(current_col > cols || current_row > rows)
					break()

				if(m1[current_row, current_col] == 0) {
					m1[current_row, current_col] = 1
					current_row = current_row + 1
				}
				
				current_col = current_col + 1
			}
		}
	}
	
	if(direction == "backward" || direction == "both") {
		rows = sum(parent)
		if(rows > 1) { # skipped if only 0s
			m2 = matrix(parent, ncol = cols, nrow = rows, byrow = TRUE)

			current_row = 1
			current_col = 1
			repeat {
				if(current_col > cols || current_row > rows)
					break()

				if(m2[current_row, current_col] == 1) {
					m2[current_row, current_col] = 0
					current_row = current_row + 1
				}
				
				current_col = current_col + 1
			}
		}
	}
	
	m = rbind(m1, m2)
	if(is.null(m))
		return(m)
	if(!is.null(omit.func)) {
		rows_to_omit = apply(m, 1, omit.func)
		return(m[!rows_to_omit,, drop = FALSE])
	} else {
		return(m)
	}
}

find.best <- function(results, subset = rep(TRUE, length(results))) {
	best = list(
		result = NULL,
		idx = NULL
	)

	w = which(subset)
	if(length(w) > 0) {
		children_results = results[w]
		max_idx = which.max(children_results)
		best$result = children_results[max_idx]
		best$idx = w[max_idx]
	}
	return(best)
}
