simulateT <- function(smp,n,N,grd,R)
{
  
  if (length(smp) != length(n) )  stop("smp and n not the same length")
  if (length(n) != length(N) )  stop("n and N not the same length")
  if ( any(round(smp) != smp) || min(smp) < 0) 
	stop("data are not  positive integers")
  if ( any(round(n) != n) || min(n) < 1) 
	stop("sample sizes are not  positive integers")
  if ( any(round(N) != N) || min(N) < 1) 
	stop("strata sizes are  not  positive integers")
  if ( min(n-smp) < 1 ) stop("data smp are  too big")
  if ( min(N-n) < 1 ) stop("sample sizes are  too big")



	K <- length(smp)
	L <- length(grd)
	dum <- rep(0,L)
	ans <- rep(0,R)

	for(k in 1:K){
		for(i in 1:L){
			dm <- gamma(smp[k] + grd[i])*gamma(n[k]-smp[k] + 1-grd[i])
			dum[i] <- dm/(gamma(grd[i])*gamma(1-grd[i]))
		}

		wprob <- dum/sum(dum)
		cumwprob <- cumsum(wprob)

			for(j in 1:R){
				v <- runif(1)
				m <- length(cumwprob[cumwprob < v])+1
				grdv <- grd[m] 
				pp <- rbeta(1,smp[k] + grdv, n[k]-smp[k] + 1-grdv)
				ans[j] <- ans[j] + smp[k] + rbinom(1,N[k]-n[k],pp)
			} 
	}

	return(ans)

}



