covI <-
function(II,m,n,ll, ThePsiJ){

K <- length(II)

Pll <- ThePsiJ[[ll]]
Nll <- (length(Pll)-1)/2

bigans <- 0

for(k in 1:K)	{

	Pk <- ThePsiJ[[k]]
	Nk <- (length(Pk)-1)/2

	mintau <- max(-Nll+n-m, -Nk)
	maxtau <- min(Nll+n-m, Nk)

	#cat("mintau:maxtau", mintau, maxtau, "\n")

	if (mintau <= maxtau)	{
		v <- mintau:maxtau
		ans <- sum(Pk[v+Nk+1]*Pll[m-n+v+Nll+1])
		}
	else
		ans <- 0

	#cat(k, ans, "\n")

	bigans <- bigans + II[k]*ans
	}
2*bigans^2
}
