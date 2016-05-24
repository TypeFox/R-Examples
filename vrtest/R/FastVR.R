FastVR <-
function(x,kvec)
{
x <- as.matrix(x)
m <- acf(x,lag.max = max(kvec),plot=FALSE)$acf[2:(max(kvec)+1)]  
VR <- matrix(NA,nrow=length(kvec))
for(i in 1:length(kvec))
{
k <- kvec[i]
w <- (1-(1:(k-1))/k)
VR[i] <- 1+ 2*sum(w*m[1:(k-1)])
}
return(VR)
}
