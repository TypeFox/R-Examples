newjk <- function(x, y, n, jlo, jhi, klo, khi, sigma2, tau2)
{
	c <- (-10^120)
	fjk <- matrix(c, n, n)
	
	##Yulei's edit to avoid using for-loop
	jkgrid <- expand.grid(jlo:jhi, klo:khi)
	res <- data.frame(j = jkgrid[,1],
	                  k = jkgrid[,2],
	                  k.ll = apply(jkgrid, 1, lpjk, x = x, 
	                               y = y, n = n,sigma2=sigma2,tau2=tau2))
	
	res.m <- matrix(res$k.ll, nrow = jhi-jlo+1, ncol = khi-klo+1)
	rownames(res.m) <- jlo:jhi
	colnames(res.m) <- klo:khi
	fjk[jlo:jhi,klo:khi]<-res.m

	z <- findmax(fjk)
	jcrit <- z$imax
	kcrit <- z$jmax
	list(jhat = jcrit, khat = kcrit, value = max(fjk))
}