transLS <- function(object, s, t, h, nh=40, ncv=10, window="normal", state.names=c("1", "2", "3"), conf=FALSE, n.boot=1000, conf.level=0.95, method.boot="percentile", boot.cv=FALSE, cv.full=TRUE) {
	if ( missing(object) ) stop("Argument 'object' is missing, with no default")
	if ( with( object[[1]], any(Stime == time1 & event) ) ) warning("Your data suggests an illness-death model\n  LS estimator is appropriate for the three-state progressive model only")
	if ( missing(s) ) s <- 0
	if ( missing(t) ) t <- max(object[[1]]$Stime)
	Message <- TPWindowStateBootCvalCheck(object, s, t, h, nh, ncv, window, state.names, conf, n.boot, conf.level, method.boot, boot.cv, cv.full)
	if ( !is.null(Message) ) stop(Message)
	if (length(h) == 1) h <- rep(h, 4)
	else if (length(h) == 2) h <- rep(sort(h), 2)
	else if (length(h) == 3) h <- c( rep(h[1], 2), sort(h[2:3]) )
	else h <- c( sort(h[1:2]), sort(h[3:4]) )
	TransLS(object)
	return( TransMethod(object, s, t, state.names, conf, n.boot, conf.level, method.boot, h=h, nh=nh, ncv=ncv, window=window, bootcv=boot.cv, cvfull=cv.full) )
}
