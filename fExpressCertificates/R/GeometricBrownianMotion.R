# Simulates one trajectory of the Geometric Brownian Motion
# at time intervals dt=T/N
# 
# @param S0 start value at t0
# @param mu Erwartungswert
# @param sigma volatility
# @param T time in years
# @return vektor of length (N+1) with prices at (i*T/N) for i=0..N.
GeometricBrownianMotion <- function(S0, mu, sigma, T, N) {
	# time intervall, grid size
	dt <- T/N
	
	S_t <- numeric(N+1)
	S_t[1] <- S0
	
	# Wiener Prozess
	W_t  <- rnorm(N,0,1)
	
	S_t[2:(N+1)] <- S0 * exp(cumsum((mu - sigma^2/2)*dt  + sigma * sqrt(dt)*W_t))
	S_t
}
# alias for method
GBM <- GeometricBrownianMotion

# Simulates multiples trajectories of the Geometric Brownian Motion.
# Returns a (mc.loops x N+1) matrix
# Ist das schneller als eine Schleife um GBB? Nein, nicht wirklich!
#
# 
# @param S0 Startwert Aktie
# @param mu Erwartungswert
# @param sigma
# @param T
# @param mc.loops number of trajectories to sample
# @return Matrix (mc.loops x N+1)
GeometricBrownianMotionMatrix <- function(S0, mu, sigma, T, mc.loops, N)
{
	# Zeitintervall
	dt <- T/N
	
	S_t     <- matrix(mc.loops, N+1)
	S_t[,1] <- S0
	
	# Wiener Prozess als Matrix (mc.loops x N)
	W_t  <- matrix(rnorm(mc.loops * N,0,1), mc.loops, N)
	
	X <- (mu - sigma^2/2)*dt  + sigma * sqrt(dt)*W_t
	
	S_t[,2:(N+1)] <- S0 * t(exp(apply(X, 1, cumsum)))
	S_t
}
GBMMatrix <- GeometricBrownianMotionMatrix