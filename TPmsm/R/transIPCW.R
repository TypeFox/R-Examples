transIPCW <- function(object, s, t, x, bw="dpik", window="normal", method.weights="NW", state.names=c("1", "2", "3"), conf=FALSE, n.boot=1000, conf.level=0.95, method.boot="percentile", method.est=1, ...) {
	if ( missing(object) ) stop("Argument 'object' is missing, with no default")
	if ( missing(s) ) s <- 0
	if ( missing(t) ) t <- max(object[[1]]$Stime)
	if (missing(x) | ncol(object[[1]]) < 5) Message <- TPStateBootCheck(object, s, t, state.names, conf, n.boot, conf.level, method.boot)
	else Message <- TPCWindowStateBootCheck(object, s, t, x, bw, window, method.weights, state.names, conf, n.boot, conf.level, method.boot)
	if ( !is.null(Message) ) stop(Message)
	if ( !is.numeric(method.est) ) stop("Argument 'method.est' must be numeric")
	if ( !(method.est %in% 1:2) ) stop("Argument 'method.est' must be one of 1 or 2")
	if (missing(x) | ncol(object[[1]]) < 5) {
		x <- NULL
		TransIPCW1(object)
	} else TransIPCW2(object)
	return( TransMethod(object, s, t, state.names, conf, n.boot, conf.level, method.boot, xi=x, bw=bw, window=window, methodweights=method.weights, methodest=method.est, ...) )
}
