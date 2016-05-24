"qvaluebh95" <-
function(p,fdrate=0.1) {
	ps <- sort.int(p,index.return=T)
	ll <- length(p)
	mul <- fdrate*c(1:ll)/ll
	pass <- (ps$x <= mul)
	if (length(ps$x[pass])) pmin <- max(ps$x[pass]) else pmin <- -1
	out <- list()
	out$significant <- (p <= pmin)
	aaa <- .C("comp_qval",as.double(ps$x),as.integer(ll), out = double(ll), PACKAGE = "GenABEL")$out
	out$qvalue[ps$ix] <- aaa
	out
}

