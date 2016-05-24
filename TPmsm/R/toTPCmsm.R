toTPCmsm <- function(lst, UT, UX, s, t, x, statenames) {
	return( .Call(Rf_toTPCmsm, lst, UT, UX, s, t, x, statenames, PACKAGE="TPmsm") )
}
