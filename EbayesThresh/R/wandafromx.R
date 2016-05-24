"wandafromx" <-
function(x)
{
#  finds the marginal max lik estimators of w and a, using a bivariate optimization
#
#  The threshold is constrained to lie between 0 and sqrt ( 2 log (n))
#
#  If running R, the routine optim is used; in S-PLUS the routine is nlminb
#
	thi <- sqrt(2 * log(length(x)))
	lo  <-  c(0,0.04)
	hi  <-  c(thi,3)
	startpar  <-  c(1,0.5)
	if (exists("optim")) {
 		uu <- optim(startpar, negloglik.laplace, method="L-BFGS-B",
			lower = lo, upper = hi, xx = x)
               	uu <- uu$par
		}
	else {uu <- nlminb(startpar, negloglik.laplace, lower = lo, upper = hi, xx = x)
	uu <- uu$parameters}
	a <- uu[2]
	w <- wfromt(uu[1], a = a)
	return(list(w=w, a=a))
}
