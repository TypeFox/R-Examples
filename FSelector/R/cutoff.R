cutoff.k <- function(attrs, k) {
	if(dim(attrs)[1] == 0)
		return(character(0))
	if(k < 1)
		stop("k too small")
	if(k > dim(attrs)[1])
		k = dim(attrs)[1]
	sorted_names = rownames(attrs)[order(attrs, decreasing = TRUE)]
	return(sorted_names[1:k])
}

cutoff.k.percent <- function(attrs, k) {
	if(dim(attrs)[1] == 0)
		return(character(0))
	if(k <= 0)
		stop("k too small")
	if(k > 1) {
		warning("Assumed k=1")
		k = 1
	}
	sorted_names = rownames(attrs)[order(attrs, decreasing = TRUE)]
	return(sorted_names[1:round((k * length(sorted_names)))])
}

cutoff.biggest.diff <- function(attrs) {
	if(dim(attrs)[1] == 0)
		return(character(0))
	else if(dim(attrs)[1] == 1)
		return(dimnames(attrs)[[1]])
	
	perm = order(attrs[,1], decreasing = TRUE)
	attrs = attrs[perm, , drop = FALSE]
	
	intervals = sapply(1:(dim(attrs)[1] - 1), function(idx) {
			attrs[idx, 1] - attrs[idx + 1, 1]
		})
	
	return(dimnames(attrs)[[1]][1:(which.max(intervals))])
}
