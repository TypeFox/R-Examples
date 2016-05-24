newst <- function(x, y, n, j0, k0, rr, ss, vv, ww)
{
	g <- rss(x, y, n, j0, k0)
	rss1 <- g$S1
	rss2 <- g$S2
	rss3 <- g$S3
	s2 <- (2 * ss + rss1 + rss3)/(n - k0 + j0 + 2 * rr + 2)
	t2 <- (2 * ww + rss2)/(k0 - j0 + 2 * vv + 2)
	list(sig2hat = s2, tau2hat = t2)
}
