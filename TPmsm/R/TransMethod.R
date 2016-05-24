TransMethod.TP <- function(object, s, t, state.names, conf, n.boot, conf.level, method.boot, ...) {
	if (conf) return( TransBoot(object, s, t, state.names, n.boot, conf.level, method.boot, ...) )
	UT <- uniqueTIME(object, s, t)
	lst <- TransPROB(object, UT, nboot=1, ...)
	class(lst) <- class(object)[2]
	TPmsm <- toTPmsm(lst, UT, s, t, state.names)
	Clean(object)
	return(TPmsm)
}

TransMethod.TPC <- function(object, s, t, state.names, conf, n.boot, conf.level, method.boot, xi, ...) {
	if (conf) return( TransBoot(object, s, t, state.names, n.boot, conf.level, method.boot, xi, ...) )
	UT <- uniqueTIME(object, s, t)
	UX <- uniqueCOV(object, xi)
	lst <- TransPROB(object, UT, UX, nboot=1, ...)
	class(lst) <- class(object)[2]
	TPCmsm <- toTPCmsm(lst, UT, UX, s, t, xi, state.names)
	Clean(object)
	return(TPCmsm)
}
