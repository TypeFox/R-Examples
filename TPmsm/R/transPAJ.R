transPAJ <- function(object, s, t, state.names=c("1", "2", "3"), conf=FALSE, n.boot=1000, conf.level=0.95, method.boot="percentile") {
	if ( missing(object) ) stop("Argument 'object' is missing, with no default")
	if ( missing(s) ) s <- 0
	if ( missing(t) ) t <- max(object[[1]]$Stime)
	Message <- TPStateBootCheck(object, s, t, state.names, conf, n.boot, conf.level, method.boot)
	if ( !is.null(Message) ) stop(Message)
	TransPAJ(object)
	return( TransMethod(object, s, t, state.names, conf, n.boot, conf.level, method.boot) )
}
