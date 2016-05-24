"dhhrcomp" <-
function (nsims=1024,nmu=4,mu0=c(0,5,10,15,20,25,30,40,50,60,65,70,80,100,200,300,500,1000)*1000) 
{
mu0 <- mu0[1:nmu]
v.hft <- rep(0,nmu)	# Store for variance for hft method
s.hft <- rep(0,nmu)	# Store for skewness
k.hft <- rep(0,nmu)	# Store for kurtosis

# Setup DHHR parameters
#
alpha <- 24800
seta <- 0.227
seps <- 4800
Seta <- sqrt(exp(seta^2)*( exp(seta^2)-1))
cc <- seps^2 / Seta^2

mu <- rep(mu0, rep(nsims, nmu))

y <- simdurbin2(mu=mu,alpha=alpha,seta=seta, seps=seps)
ddhfty <- ddhft.np.2(y)

for(i in 1:nmu)	{
	sv <- mu==mu0[i]
	v.hft[i] <- var(ddhfty$hft[sv])
	s.hft[i] <- SKEW(ddhfty$hft[sv])
        k.hft[i] <- KURTOSIS(ddhfty$hft[sv])
	}
	
l <- list(mu=mu0, v.hft=v.hft, s.hft=s.hft, k.hft=k.hft)
return(l)


}

