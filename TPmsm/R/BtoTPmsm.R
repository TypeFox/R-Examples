BtoTPmsm.default <- function(lst, UT, s, t, statenames, nboot, conflevel, methodboot) {
	return( .Call(Rf_BtoTPmsm1222, lst, UT, s, t, statenames, as.integer(nboot), conflevel, methodboot, PACKAGE="TPmsm") )
}

BtoTPmsm.KMW2 <- function(lst, UT, s, t, statenames, nboot, conflevel, methodboot) {
	return( .Call(Rf_BtoTPmsm1323, lst, UT, s, t, statenames, as.integer(nboot), conflevel, methodboot, PACKAGE="TPmsm") )
}

BtoTPmsm.KMPW2 <- function(lst, UT, s, t, statenames, nboot, conflevel, methodboot) {
	return( .Call(Rf_BtoTPmsm1323, lst, UT, s, t, statenames, as.integer(nboot), conflevel, methodboot, PACKAGE="TPmsm") )
}
