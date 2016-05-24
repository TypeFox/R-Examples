transKMW <- function(object, s, t, state.names=c("1", "2", "3"), conf=FALSE, n.boot=1000, conf.level=0.95, method.boot="percentile", method.est=3) {
	if ( missing(object) ) stop("Argument 'object' is missing, with no default")
	if ( missing(s) ) s <- 0
	if ( missing(t) ) t <- max(object[[1]]$Stime)
	Message <- TPStateBootCheck(object, s, t, state.names, conf, n.boot, conf.level, method.boot)
	if ( !is.null(Message) ) stop(Message)
	if ( !is.numeric(method.est) ) stop("Argument 'method.est' must be numeric")
	if ( !(method.est %in% 1:4) ) stop("Argument 'method.est' must be one of 1, 2, 3 or 4")
	if ( method.est %in% c(1, 2) ) TransKMW1(object)
	else TransKMW2(object)
	return( TransMethod(object, s, t, state.names, conf, n.boot, conf.level, method.boot, methodest=method.est) )
}
