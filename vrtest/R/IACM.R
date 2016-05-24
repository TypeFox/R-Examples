IACM <-
function(r,uhat)
{
T <- length(uhat)
iacm <- numeric(3)
seq <- 0:(T-1)

gamz <- sum(uhat^2)/T

fhat <- 0
if (r > 0)
{

for (s in 1:as.integer(T*r/2))
{    
    lamda <- 2*pi*s/T
    co <- cos(lamda*seq)
    si <- sin(lamda*seq)
    perio <- ( (sum(co*uhat) )^2 + ( sum(si*uhat) )^2 )/(2*pi*T)
    fhat <- fhat + perio
}

}

fhat <- 2*pi*fhat/T

ut <- fhat - gamz*r/2
uts <- sqrt(2*T)*ut/gamz

iacm[1] <- uts^2/(r*(1-r))
if (r == 0 | r ==1 )
iacm[1] <- 0

iacm[2] <- uts^2
iacm[3] <- abs(uts)
return(iacm)
}
