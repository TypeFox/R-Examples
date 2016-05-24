# Simulates a Brownian Bridge B(t) for 0<=t<=T with conditions B(t)=a and B(T)=b
#
BrownianBridge<-function(a=0, b=0, sigma=1, t0=0, T=1, N)
{
  # Zeitpunkte
  t <- seq(t0, T, by=(T-t0)/N)
  
  # Brownsche Bewegung ohne Drift
  W_t <- BrownianMotion(S0=a, mu=0, sigma=sigma, T=T, N=N)
  
  # Brownsche Brücke
  B_t <- W_t - ((t-t0)/(T-t0))*W_t[N+1] + ((t-t0)/(T-t0)) * b
  
  B_t
}

# Draws from the minimum distribution m(T) of a Brownian Bridge B(t) für 0 <= t <= T mit T>0, B(0)=0, B(T)=a 
# (see Beskos et al. (2004), p.7,  generalized to the sigma^2<>1 case)
#
# Simulation uses the conditional density f(m(T) | B(T)=a) (see Beskos et al. (2004), p.20), 
# scaled by sigma^2.
#
# @param n number of samples to create
# @param t0
# @param T
# @param a B(t0) = a
# @param b B(T) = b
# @return a vector of length n with samples from the minimum
rBrownianBridgeMinimum <- function(n=100, t0=0, T=1, a=0, b=0, sigma=1)
{
  # Exponential verteilte Zufallsgröße
  E <- rexp(n)
  
  # Verschiebung von W(t_0)=a und W(T)=b auf W(t_0)=0 und W(T)=b-a
  a_norm <- b-a
  
  # Skalierung in der Zeitachse
  T <- sigma^2 * (T-t0)
  
  # Simuliere das Minimum m_T | W_T = a
  Z1 <- (a_norm-sqrt(2*T*E + a_norm^2))/2
  
  # Rücktransformation (Verschiebung um +a)
  Z1+a
}

# density function of the minimum m_T of a Brownian Bridge B_t
# for 0 <= t <= T and T>0, B_{t0}=a, B_T=b,  
# which is the conditional density f(m_T = x | B_{t0}=a, B_T=b)
# (see Beskos et al. (2004), p.7/p.20,  generalized to the sigma^2<>1 case)
#
# @param t0
# @param T
# @param a B(t0) = a
# @param b B(T) = b
#
dBrownianBridgeMinimum <- function(x, t0=0, T=1, a=0, b=0, sigma=1)
{
  t <- (T-t0)
  2 * (a-2*x) / (sigma^2 * t) * exp(- ((a - 2*x)^2 - a^2)/(2*sigma^2*t))
}
  
