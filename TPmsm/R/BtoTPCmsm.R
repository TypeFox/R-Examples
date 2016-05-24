BtoTPCmsm <- function(lst, UT, UX, s, t, x, statenames, nboot, conflevel, methodboot) {
	return( .Call(Rf_BtoTPCmsm, lst, UT, UX, s, t, x, statenames, as.integer(nboot), conflevel, methodboot, PACKAGE="TPmsm") )
}
