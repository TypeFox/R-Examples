TPBoot.rboot <- function(object, UT, nboot, ...) {
	lst <- vector(mode="list", length=2)
	a3d <- array( data=0, dim=c(nboot, length(UT), 4) )
	est <- TransPROB(object, UT, nboot=1, ...)
	a3d[1,,] <- est[[1]]
	lst[[2]] <- est[[2]]
	for (i in 2:nboot) {
		y <- vector(mode="list", length=1)
		y[[1]] <- object[[1]][csample.int( n=nrow(object[[1]]) ),] # csample.int must be called for reproduction with TPmsm's own RNG
		class(y) <- class(object)[2]
		a3d[i,,] <- TransPROB(y, UT, nboot=1, ...)[[1]]
	}
	lst[[1]] <- a3d
	return(lst)
}

TPBoot.cboot <- function(object, UT, nboot, ...) {
	return( TransPROB(object, UT, nboot=nboot, ...) )
}

TransBoot.TP <- function(object, s, t, state.names, n.boot, conf.level, method.boot, ...) {
	UT <- uniqueTIME(object, s, t)
	n.boot <- trunc(n.boot)
	lst <- TPBoot(object, UT, n.boot, ...)
	class(lst) <- class(object)[2]
	TPmsm <- BtoTPmsm(lst, UT, s, t, state.names, n.boot, conf.level, method.boot)
	Clean(object)
	return(TPmsm)
}

TPCBoot.rboot <- function(object, UT, UX, nboot, ...) {
	lst <- vector(mode="list", length=2)
	a4d <- array( data=0, dim=c(nboot, length(UT), length(UX), 4) )
	est <- TransPROB(object, UT, UX, nboot=1, ...)
	a4d[1,,,] <- est[[1]]
	lst[[2]] <- est[[2]]
	for (i in 2:nboot) {
		y <- vector(mode="list", length=1)
		y[[1]] <- object[[1]][csample.int( n=nrow(object[[1]]) ),] # csample.int must be called for reproduction with TPmsm's own RNG
		class(y) <- class(object)[2]
		a4d[i,,,] <- TransPROB(y, UT, UX, nboot=1, ...)[[1]]
	}
	lst[[1]] <- a4d
	return(lst)
}

TPCBoot.cboot <- function(object, UT, UX, nboot, ...) {
	return( TransPROB(object, UT, UX, nboot=nboot, ...) )
}

TransBoot.TPC <- function(object, s, t, state.names, n.boot, conf.level, method.boot, xi, ...) {
	UT <- uniqueTIME(object, s, t)
	UX <- uniqueCOV(object, xi)
	n.boot <- trunc(n.boot)
	lst <- TPCBoot(object, UT, UX, n.boot, ...)
	class(lst) <- class(object)[2]
	TPCmsm <- BtoTPCmsm(lst, UT, UX, s, t, xi, state.names, n.boot, conf.level, method.boot)
	Clean(object)
	return(TPCmsm)
}
