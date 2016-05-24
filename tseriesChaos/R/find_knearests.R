#Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-02 17:15:39 +0100 (ven, 02 dic 2005) $
find_knearests <- function(series, m, d, t, eps, ref, k, s){
	x <- ( series - min(series) ) / (diff(range(series)))
	eps <- eps/diff(range(series))
	res <- numeric(ref*k)
	res <- .C("find_knearests", series=as.double(x), m=as.integer(m), d=as.integer(d), t=as.integer(t), length=as.integer(length(x)), eps=as.double(eps),ref=as.integer(ref), k=as.integer(k), s=as.integer(s), res=as.integer(res), PACKAGE="tseriesChaos")[["res"]]
	res[res==-1] <- NA;
	matrix(res, ref, k)
}
