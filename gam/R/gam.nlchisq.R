"gam.nlchisq" <-
function(qr, resid, w, s)
{
	wt <- sqrt(w)
	s <- s * wt
	resid <- wt * resid
	Rsw <- qr.resid(qr, s)
	apply(Rsw^2 + 2 * s * resid, 2, sum)
}
