sqndwd <-
function (x, ec) 
{
    type <- "station"
    m0 <- ec$m0
    filter.number <- ec$filter.number
    family <- ec$family
    ordwd <- wd(x, filter.number = filter.number, family = family, 
        type = type)
    answd <- wd(rep(0, length(x)), filter.number = filter.number, 
        family = family, type = type)
    nlev <- nlevelsWT(ordwd)
    if (m0 < nlev)	{
	    for (j in 0:(nlev - m0 - 1)) {
		cjm0 <- accessC(ordwd, level = j + m0)
		avec <- rep(0, length = lc <- length(cjm0))
		for (k in 0:(lc - 1)) {
		    loclvec <- k + (ec$ll) * 2^(nlev - (j + m0))
		    cix <- loclvec%%lc
		    avec[k + 1] <- sum(ec$ecoef * cjm0[cix + 1])
		}
		answd <- putD(answd, level = j, v = avec)
	    }
	}
    brute <- sqndwdecomp(x, J = m0, filter.number = filter.number, 
        family = family)
    for (j in (nlev - m0):(nlev - 1)) answd <- putD(answd, level = j, 
        v = brute[nlev - j, ])
    answd
}
