reverse.kronecker <-
function(ab,p,qq) {
	m <- nrow(ab)/p
	n <- ncol(ab)/qq
	rr <- c(outer((0:(p-1))*m,1:m,"+"))
	cc <- c(outer((0:(qq-1))*n,1:n,"+"))
	ab[rr,cc]
}
