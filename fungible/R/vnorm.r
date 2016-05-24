# scale vector to unit length
	vnorm <- function(x) as.matrix(x/as.numeric(sqrt(t(x) %*% x)))