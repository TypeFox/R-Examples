
# Simulates one trajectory of the Arithmetic Brownian Motion
# at grid points with time interval dt=T/N from (0,...,T)
#
# @param S0 Startwert Aktie
# @param mu Erwartungswert
# @param sigma
# @param T Laufzeit
# @param N number of grid points
#
# @return S a grid of (N+1) samples including start value S0, i.e. N simulates samples
BrownianMotion<-function(S0, mu=0, sigma=1, T, N)
{
	# time interval
	dt <- T/N
	
	# Speicher reservieren und S_0 belegen
	S_t <- numeric(N+1)
	S_t[1] <- S0
	
	# Wiener Prozess
	W_t <- rnorm(N,0,1)
	
	S_t[2:(N+1)] <- S0 + cumsum(mu*dt  + sigma * sqrt(dt)*W_t)
	S_t
}

# alias
BM <- BrownianMotion





