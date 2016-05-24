`toParamIndexes` <-
function(m, string) {
	if (is.character(m)) {
		ans= matchIndexes(m, string)
		if (any(is.na(ans))) {
			stop("error: at least one parameter name cannot be found")
		} 
		return (ans)	
	}
	return (m)
}

