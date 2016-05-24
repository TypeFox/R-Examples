"PoincarePlot" <-
function(x, res)
{
  ind <- order(x)
  e <- res[ind]
	et <- e[ - length(e)]
	etp1 <- e[-1]
	plot(et, etp1, xlab = "e[t]", ylab = "e[t+1]")
	lines(lowess(et, etp1, f = 1), lwd = 2)
	invisible()
}
