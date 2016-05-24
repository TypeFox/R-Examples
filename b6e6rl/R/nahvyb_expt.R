nahvyb_expt <-
function(N,k,expt){
opora <- 1:N
nargin <- length(as.list(match.call())) -1
if (nargin==3) 
	opora <- opora[-expt]
vyb <- rep(0,k)
for (i in 1:k){
	index <- 1+trunc(runif(1)*length(opora))
	vyb[i] <- opora[index]
	opora <- opora[-index]
}
return(vyb)
}
