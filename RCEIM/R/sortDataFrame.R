sortDataFrame <- function(x, key, ...) {
# Sorting a dataframe by a given key
# From: http://snippets.dzone.com/posts/show/470
# Post by: r-fanatic 
	if (missing(key)) {
		rn <- rownames(x)
		if (all(rn %in% 1:nrow(x))) rn <- as.numeric(rn)
		x[order(rn, ...), , drop=FALSE]
	} else {
		x[do.call("order", c(x[key], ...)), , drop=FALSE]
	}
}