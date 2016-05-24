llsearch <- function(x, y, n, jlo, jhi, klo, khi,plot)
{
  
	fjk <- matrix(0, n, n)
	fxy <- matrix(0, (jhi - jlo + 1), (khi - klo + 1))
	
  ## Yulei's edit to avoid using for-loop
	jkgrid <- expand.grid(jlo:jhi, klo:khi)
	k.ll <- apply(jkgrid, 1, p.estFUN, x = x, 
	                               y = y, n = n)
	
	fxy <- matrix(k.ll, nrow = jhi-jlo+1, ncol = khi-klo+1)
	rownames(fxy) <- jlo:jhi
	colnames(fxy) <- klo:khi
	if (plot == "TRUE") {
	 jx<-jlo:jhi
   ky<-klo:khi
    persp(jx, ky, fxy, xlab = "j", ylab = "k", zlab = "LL(x,y,j,k)")
	  title("Log-likelihood Surface")
	}
	z <- findmax(fxy)
	jcrit <- z$imax + jlo - 1
	kcrit <- z$jmax + klo - 1
	list(jhat = jcrit, khat = kcrit, value = max(fxy))
}
