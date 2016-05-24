TransPROB.AJ <- function(object, UT, nboot, ...) {
	return( .Call(Rf_TransPROBAJ, object, UT, as.integer(nboot), PACKAGE="TPmsm") )
}

TransPROB.PAJ <- function(object, UT, nboot, ...) {
	return( .Call(Rf_TransPROBPAJ, object, UT, as.integer(nboot), PACKAGE="TPmsm") )
}

TransPROB.KMW1 <- function(object, UT, nboot, methodest, ...) {
	return( .Call(Rf_TransPROBKMW, object, UT, as.integer(nboot), as.integer(methodest), PACKAGE="TPmsm") )
}

TransPROB.KMW2 <- function(object, UT, nboot, methodest, ...) {
	return( .Call(Rf_TransPROBKMW, object, UT, as.integer(nboot), as.integer(methodest), PACKAGE="TPmsm") )
}

TransPROB.KMPW1 <- function(object, UT, nboot, methodest, ...) {
	return( .Call(Rf_TransPROBKMPW1, object, UT, as.integer(nboot), as.integer(methodest), PACKAGE="TPmsm") )
}

TransPROB.KMPW2 <- function(object, UT, nboot, methodest, ...) {
	return( .Call(Rf_TransPROBKMPW2, object, UT, as.integer(nboot), as.integer(methodest), PACKAGE="TPmsm") )
}

TransPROB.IPCW1 <- function(object, UT, nboot, methodest, ...) {
	return( .Call(Rf_TransPROBIPCW1, object, UT, as.integer(nboot), as.integer(methodest), PACKAGE="TPmsm") )
}

TransPROB.IPCW2 <- function(object, UT, UX, bw, window, methodweights, nboot, methodest, ...) {
	h <- TransWidth(object, bw, window, ...)
	return( .Call(Rf_TransPROBIPCW2, object, UT, UX, h, window, methodweights, as.integer(nboot), as.integer(methodest), PACKAGE="TPmsm") )
}

TransPROB.LIN1 <- function(object, UT, nboot, ...) {
	return( .Call(Rf_TransPROBLIN1, object, UT, as.integer(nboot), PACKAGE="TPmsm") )
}

TransPROB.LIN2 <- function(object, UT, UX, bw, window, methodweights, nboot, ...) {
	h <- TransWidth(object, bw, window, ...)
	return( .Call(Rf_TransPROBLIN2, object, UT, UX, h, window, methodweights, as.integer(nboot), PACKAGE="TPmsm") )
}

TransPROB.LS <- function(object, UT, h, nh, ncv, window, nboot, bootcv, cvfull, ...) {
	return( .Call(Rf_TransPROBLS, object, UT, h, as.integer(nh), as.integer(ncv), window, as.integer(nboot), as.logical(bootcv), as.logical(cvfull), PACKAGE="TPmsm") )
}
