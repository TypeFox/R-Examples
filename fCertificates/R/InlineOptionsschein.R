# Bewertung von Inline-Optionsscheinen / Range of the Geometric Brownian Motion / Double-Barrier Binary Option
#
# Range a < x < b
#
# P(a < m <= M < b, x) --> fixer Payoff
#
# Literatur:
#
# (1) Sutrick, Teall, Tucker, Wei (1997): The Range of Brownian Motion Processes, The Journal of Financial Engineering, Vol.6, 31-46
# (2) Hui (1996) : One-Touch Barrier Binary Option Values, Applied Financial Economics Vol. 6, 343-346
# (3) Haug (2007), p. 180/181

#source("C:/Projects/Java/certificateserver/docs/modelle/GBB.R")
#source("F:/Promotion/Sweave/Expressstudie/GBB.R")

# Simuliere N Geometrische Brownsche Bewegungen
#S0    = 100
#mu    = 0
#sigma = 0.3
#T     = 1
#
#simProbInRange <- function(S0, a, b, mu, sigma, T, mc.loops, mc.steps)
#{
#  S     = matrix(NA, mc.loops, mc.steps+1)
#  for (i in 1:mc.loops)
#  {
#    S[i,] = GBB(S0=S0, mu, sigma, T, N=mc.steps)
#  }#
#
#  # Minimum m_t
#  m_t = apply(S, 1, min)
#  M_t = apply(S, 1, max)

#  # Wahrscheinlichkeit aus Monte Carlo, zwischen a und b zu bleiben
#  p = mean(m_t >= a & M_t <= b)
#  p
#}

# Arithmetische Brownsche Bewegung ohne Drift: Wahrscheinlichkeit, in der Range [a,b] zu bleiben
#
# P(a < m <= M < b, x)
#
probInRange <- function(x, a, b, T, sigma)
{
  # k = -Inf ... +Inf
  p = 0
  for (k in -100:100)
  {
    p = p + (pnorm((b + 2 * k * (b - a)) / (sigma*sqrt(T)))        - pnorm((a + 2 * k * (b - a)) / (sigma*sqrt(T)))) - 
            (pnorm((2*b - a + 2 * k * (b - a)) / (sigma*sqrt(T)))  - pnorm((2 * b - b + 2 * k * (b - a)) / (sigma*sqrt(T))))
  }
  return (p)
}

# Wahrscheinlichkeit ausrechnen, zwischen 80 und 120 zu bleiben
#probInRange(a=log(70)-log(S0), b=log(120)-log(S0), T=1, sigma=0.3)
#
#probInRange(a=log(4800)-log(5473), b=log(5600)-log(5473), T=7/365, sigma=0.3)
#
#simProbInRange(S0=5473, a=4800, b=5600, mu=0.01, sigma=0.28, T=7/365, mc.loops=1000, mc.steps=1000)

# Formel von Hui (1996)
#
# @param S
# @param L lower barrier
# @param U upper barrier
# @param K fix payoff
DoubleBarrierBinaryCall <- function(S, K, L, U, T, r, r_d, sigma, ratio=1, nmax=20)
{
  # i = 1 ... +Inf
  c = 0
  
  # cost-of-carry
  b     = r - r_d
  Z     = log(U/L)
  alpha = -0.5  * (2*b/sigma^2 - 1)
  beta  = -0.25 * (2*b/sigma^2 - 1) - 2*r/sigma^2
  
  for (i in 1:nmax)
  {
    c = c + 2*pi*i*K/Z^2 * ((S/L)^alpha - (-1)^i * (S/U)^alpha) / (alpha^2 + (i*pi/Z)^2) * sin(i*pi/Z*log(S/L)) * exp(-0.5*((i*pi/Z)^2 - beta)*sigma^2*T)
  }
  return(c*ratio)
}

# Beispiel: CM2KC7
DoubleBarrierBinaryCall(S=5473, K=10, L=4800, U=5600, T=6/365, r=0.01, r_d=0, sigma=0.28)