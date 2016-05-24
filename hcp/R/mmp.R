mmp <- function(x, y, n, jlo, jhi, klo, khi, sig, tau, rr, ss, vv, ww)
  {
	# 1.  Given sigma^2 & tau^2, find j & k
	step1 <- newjk(x, y, n, jlo, jhi, klo, khi, sig, tau)
	jhat <- step1$jhat
	khat <- step1$khat

	# 2.  Given j & k, find sigma^2 & tau^2
	step2 <- newst(x, y, n, jhat, khat, rr, ss, vv, ww)
	sighat <- step2$sig2hat
	tauhat <- step2$tau2hat
	list(jnew = jhat, knew = khat, signew = sighat, taunew = tauhat)
}