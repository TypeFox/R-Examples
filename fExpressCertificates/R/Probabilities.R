# Some formulas for analytic probabilities for the 
# Arithmetic Brownian Motion (BM) and Geometric Brownian Motion (GBM)
#
# e.g. for maximum/minumum, joint maximum and price etc.
#
# See :
#
# Chuang (1996): Joint distribution of Brownian motion and its maximum, with a 
#                generalization to correlated BM and applications to barrier options
#                \emph{Statistics & Probability Letters} (28), pp.81-90
# Poulsen :

# Simulation der Wahrscheinlichkeit, bis zum Zeitpunkt in der Range [a,b] zu bleiben,
# d.h. P(m_t >= a, M_t <= b) = P(a <= m_t <= M_t <= b)
#
# @param S0
# @param a lower bound
# @param b upper bound
# @param mu
# @param sigma
# @param T
# @param mc.loops
# @param mc.steps
simProbInRange <- function(S0, a, b, mu, sigma, T, mc.loops, mc.steps)
{
  S <- matrix(NA, mc.loops, mc.steps+1)
  for (i in 1:mc.loops) {
    S[i,] <- GBM(S0=S0, mu, sigma, T, N=mc.steps)
  }

  # Minimum m_t
  m_t <- apply(S, 1, min)
  M_t <- apply(S, 1, max)

  # Wahrscheinlichkeit aus Monte Carlo, zwischen a und b zu bleiben
  p <- mean(m_t >= a & M_t <= b)
  p
}

# Arithmetische Brownsche Bewegung ohne Drift: Wahrscheinlichkeit, in der Range [a,b] zu bleiben
#
# P(a <= m_t <= M_t <= b, x)
#
probInRange <- function(x, a, b, T, sigma)
{
	# k = -Inf ... +Inf
	p <- 0
	for (k in -100:100)
	{
		p <- p + (pnorm((b + 2 * k * (b - a)) / (sigma*sqrt(T)))        - pnorm((a + 2 * k * (b - a)) / (sigma*sqrt(T)))) - 
				(pnorm((2*b - a + 2 * k * (b - a)) / (sigma*sqrt(T)))  - pnorm((2 * b - b + 2 * k * (b - a)) / (sigma*sqrt(T))))
	}
	return (p)
}

# Calculate probabilities for the Arithmetic Brownian Motion B_t with drift parameter mu and volatility sigma
#
# @param Type different types of probabilities
# @param a
# @param z
# @param t
# @param mu drift parameter mu
# @param sigma volatility sigma
calculateProbabilityBrownianMotion <- function(Type=
c("P(M_t >= a)",
  "P(M_t <= a)",
  "P(m_t <= a)",
	"P(m_t >= a)",
  "P(M_t >= a, B_t <= z)",
  "P(m_t <= a, B_t >= z)",
  "P(a <= m_t, M_t <= b)",
  "P(M_s >= a, B_t <= z | s < t)",
  "P(m_s <= a, B_t >= z | s < t)",
  "P(M_s >= a, B_t <= z | s > t)",
  "P(m_s <= a, B_t >= z | s > t)"), a, z=0, t=1, mu=0, sigma=1, s=0)
{
  Type <- match.arg(Type)
  
  if (Type == "P(M_t >= a)") {
    # Maximum greater than a, P(M_t >= a), Chuang (1996), Gleichung (2.3)
    1 - pnorm((a - mu*t) / (sigma*sqrt(t))) + exp(2*mu*a/sigma^2) * pnorm((-a - mu*t) / (sigma*sqrt(t)))
  } else if (Type == "P(M_t <= a)") {
    # Maximum smaller than a, P(M_t <= a) = 1 - P(M_t >= a), Chuang (1996), Gleichung (2.3)
	  pnorm((a - mu*t) / (sigma*sqrt(t))) - exp(2*mu*a/sigma^2)      * pnorm((-a - mu*t) / (sigma*sqrt(t)))
  } else if (Type == "P(m_t <= a)") {
    # P(m_t <= a), Chuang (1996), Gleichung (2.3), aber für das Minimum, a < 0
    1 - pnorm((-a + mu*t) / (sigma*sqrt(t))) + exp(2*mu*a/sigma^2) * pnorm((a + mu*t) / (sigma*sqrt(t)))	
  } else if (Type == "P(m_t >= a)") {
	  # P(m_t >= a) = 1 - P(m_t <= a), Chuang (1996), Gleichung (2.3), aber für das Minimum, a < 0
	  pnorm((-a + mu*t) / (sigma*sqrt(t))) - exp(2*mu*a/sigma^2) * pnorm((a + mu*t) / (sigma*sqrt(t)))	
  } else if (Type == "P(M_t >= a, B_t <= z)") {
    # P(M_t >= a, B_t <= z) vgl. Chuang (1996), Gleichung (2.1), S.82
    if (a < 0 ) stop("a muss groesser sein als 0")
	
	  if (z <= a) {
	    exp(2*mu*a/sigma^2) * pnorm((z - 2*a - mu*t)/(sigma*sqrt(t)))
	  } else {
	  # z > a	
	    pnorm((z - mu*t)/(sigma*sqrt(t))) - pnorm((a - mu*t)/(sigma*sqrt(t))) + exp(2*mu*a/sigma^2) * pnorm((-a - mu*t)/(sigma*sqrt(t)))
	  }
  } else if (Type == "P(m_t <= a, B_t >= z)") {
    # joint.cdf.m_t.B_t
    # P(m_t <= a, B_t >= z) vgl. Chuang (1996), Gleichung (2.1), S.82, für das Minimum
    if (a > 0) stop("a muss kleiner sein als 0")
	  if (z >= a) {
	    exp(2*mu*a/sigma^2) * pnorm((-z + 2*a + mu*t)/(sigma*sqrt(t)))
	  } else {
	  # z < a	
	    pnorm((-z + mu*t)/(sigma*sqrt(t))) - pnorm((-a + mu*t)/(sigma*sqrt(t))) + exp(2*mu*a/sigma^2) * pnorm((a + mu*t)/(sigma*sqrt(t)))
	  }
  } else if (Type == "P(M_s >= a, B_t <= z | s < t)") {
    # Formel (2.7) nach Chuang (1996), S.84 für das gemeinsame Maximum M_s und den Wert B_t (M_s,B_t) s>t:
    # P(M_s >= a, B_t <= z) mit s<t
    joint.cdf.M_s.B_t1(a=a, z=z, s=s, t=t, mu=mu, sigma=sigma)
  } else if (Type == "P(m_s <= a, B_t >= z | s < t)") {
    # Formel (2.7) nach Chuang (1996), S.84 für das Minimum m_s und den Wert B_t : 
    # P(m_s <= a, B_t >= z) mit s<t    (a < 0, z < 0, s < t)
    joint.cdf.m_s.B_t1(a=a, z=z, s=s, t=t, mu=mu, sigma=sigma)
  } else if (Type == "P(M_s >= a, B_t <= z | s > t)") {
    # Formel (2.9) nach Chuang (1996), S.85 für das Maximum M_s und den Wert B_t mit s>t : 
    # P(M_s >= a, B_t <= z) mit s>t
    joint.cdf.M_s.B_t(a=a, z=z, s=s, t=t, mu=mu, sigma=sigma)
  } else if (Type == "P(m_s <= a, B_t >= z ) | s > t)") {
    # Formel (2.9) nach Chuang (1996), S.85 für die gemeinsame Wahrscheinlichkeit (M_s, B_t) des Maximum M_s und den Wert B_t mit s>t : 
    # P(M_s >= a, B_t <= z) mit s>t
    #
    # Für das Minimum gilt dann
    # --> P(m_s <= -a, B*_t >= -z)
    #
    joint.cdf.M_s.B_t(a=a, z=z, s=s, t=t, mu=mu, sigma=sigma)
  }
}
calcBMProbability <- calculateProbabilityBrownianMotion

# Formel (2.7) nach Chuang (1996), S.84 für das Maximum M_s und den Wert B_t (M_s,B_t): 
#
# P(M_s >= a, B_t <= z) mit s<t
#
joint.cdf.M_s.B_t1 <- function(a=110, z=100, s=0.6, t=1, mu=0.05, sigma=1)
{
  if (s>=t) stop("Must be s < t")
  
	# Korrelation zwischen den Zeitpunkten s und t
  rho <- sqrt(s/t)
    
  # correlation matrix between s and t
  R  <- matrix(c(1,rho,rho,1), 2, 2)
    
  # 
  x <- (-a - mu*s)/(sigma*sqrt(s))
    
	p1 <- pmvnorm(upper=c(x, (z-2*a-mu*t)/(sigma*sqrt(t))), corr = R)
	p2 <- pmvnorm(upper=c((a - mu*s)/(sigma*sqrt(s)), (z    -mu*t)/(sigma*sqrt(t))), corr = R)
	
	exp(2*mu*a/sigma^2) * p1 + pnorm(z-mu*t/(sigma*sqrt(t))) - p2
}


# Formel (2.7) nach Chuang (1996), S.84 für das Minimum m_s und den Wert B_t : 
#
# P(m_s <= a, B_t >= z) mit s<t    (a < 0, z < 0, s < t)
#
joint.cdf.m_s.B_t1 <- function(a=-1, z=-0.2, s=0.6, t=1, mu=0.05, sigma=1)
{
  if (s>=t) stop("Must be s < t")
  if (a>=0) stop("Minimum must be a < 0")
  if (z>=0) stop("Must be z < 0")
  
	# Korrelation zwischen den Zeitpunkten s und t
  rho = sqrt(s/t)
    
  # correlation matrix between s and t
  R   = matrix(c(1,rho,rho,1), 2, 2)
    
  # 
  x = (a + mu*s)/(sigma*sqrt(s))
  y = (-z +2*a + mu*t)/(sigma*sqrt(t))
  x2 = (-a + mu*s)/(sigma*sqrt(s))
  y2 = (-z + mu*t)/(sigma*sqrt(t))
    
  p1 <- pmvnorm(upper=c(x, y), corr = R)
  p2 <- pnorm(-z+mu*t/(sigma*sqrt(t)))
  p3 <- pmvnorm(upper=c(x2,y2), corr = R)
	
  exp(2*mu*a/sigma^2) * p1 + p2 - p3
}

# Formel (2.9) nach Chuang (1996), S.85 für die gemeinsame Wahrscheinlichkeit (M_s, B_t) des Maximum M_s und den Wert B_t mit s>t : 
# P(M_s >= a, B_t <= z) mit s>t
#
# Für das Minimum gilt dann
# --> P(m_s <= -a, B*_t >= -z)
#
joint.cdf.M_s.B_t<-function(a=0, z=0, s=1, t=0.6, mu=0.05, sigma=1)
{
	# P(M_t >= a, B_t <= z) vgl. Gleichung (2.1), S.82
	P1 = joint.cdf.M_t.B_t(a=a, z=z, t=t, mu=mu, sigma=sigma)
    
    if (z < a)
    {
      m1 = min(z, a)
    
	  u  = (-a + mu * s) / (sigma*sqrt(s))
      v  = (-m1 + mu * t) / (sigma*sqrt(t))
      x  = -(a + mu * s) / (sigma*sqrt(s))
      y  = -(m1 + mu * t) / (sigma*sqrt(t))
      b  = 2 * a / (sigma*sqrt(t))
	
	  # Korrelation zwischen den Zeitpunkten
	  rho = sqrt(t/s)
    
      # Korrelationsmatrix R
      # correlation matrix between s and t (s > t)
      R   = matrix(c(1,rho,rho,1), 2, 2) 
	
      # SW: Geändert von pmvnorm(upper=c(-y+b, -u), corr=R)  auf pmvnorm(upper=c(y+b, -u), corr=R) - 
	  res = P1 + pnorm(u) - pnorm(-u) + 
      pmvnorm(upper=c(y+b, -u), corr=R) - 
      pmvnorm(upper=c(v   , u) , corr=R) + 
      exp(2*mu*a/sigma^2) * (pnorm(x) - pnorm(-x) + pmvnorm(upper=c(v+b, -x), corr=R) - pmvnorm(upper=c(y, x), corr=R))
	
	}
    else # z>=a : P(M_s >= a, B_t <= z) = P(M_t >= a, B_t <= z) + P(M_s >= a) - P(M_t >= a)
    {
      P2  = cdf.M_t(a, t=s, mu, sigma) # P(M_s >= a)
      P3  = cdf.M_t(a, t=t, mu, sigma) # P(M_t >= a)
      res = P1 + P2 - P3
    }
    
    return(res)
}

# ???? Formel (2.9) nach Chuang (1996), S.85 für das Minimum m_s und den Wert B_t mit s>t : 
# P(m_s <= a, B*_t >= z) mit s>t  (a < 0)
#
#
joint.cdf.m_s.B_t <- function(a=0, z=0, s=1, t=0.6, mu=0.05, sigma=1)
{
	# P(m_t <= a, B_t >= z) vgl. Gleichung (2.1), S.82
	P1 = joint.cdf.m_t.B_t(a=a, z=z, t=t, mu=mu, sigma=sigma)
	
	m1 = max(z, a)
    
	u  = (a - mu*s) / (sigma*sqrt(s))
    v  = (m1 - mu * t) / (sigma*sqrt(t))
    x  = (a + mu*s) / (sigma*sqrt(s))
    y  = (m1 + mu * t) / (sigma*sqrt(t))
    b  = -2 * a / (sigma*sqrt(t))
	
	# Korrelation zwischen den Zeitpunkten
	rho = sqrt(t/s)
    
    # Korrelationsmatrix R
    # correlation matrix between s and t (s > t)
    R   = matrix(c(1,rho,rho,1), 2, 2) 
	
	res = P1 + pnorm(u) - pnorm(-u) + pmvnorm(upper=c(y+b, -u), corr=R) - pmvnorm(upper=c(v   , u) , corr=R) + 
        exp(2*mu*a/sigma^2) * (pnorm(x) - pnorm(-x) + pmvnorm(upper=c(v+b, -x), corr=R) - pmvnorm(upper=c(y, x), corr=R))
    return(res)
}

# Totalverlustwahrscheinlichkeit für Plain-Vanilla-Optionen : Risiko, am Ende der Laufzeit aus dem Geld zu sein. Formal : P(S_T <= X)
#
# @param type : "c" = call, "p" = put
# @param S underlying price
# @param X strike price
# @param T time to maturity 
# @param r riskfree rate
# @param r_d dividend yield
# @param sigma volatility
shortfall_risk<-function(type="c", S, X, T, r, r_d, sigma)
{ 
	mu=r-r_d
	p = pnorm((log(X/S) - (mu - sigma^2/2)*T) / sqrt(sigma^2 * T), mean=0, sd=1);
	if (type == "c") {
	  # P(S_T <= X)
	  p
	} else {
	  # P(S_T >= X)	
	  (1-p)
	}
}

# Calculate probabilities for the Geometric Brownian Motion B_t with drift mu and volatility sigma
#
#
calculateProbabilityGeometricBrownianMotion <- function(Type=c(
				 "P(S_t <= X)",
				 "P(S_t >= X)",
				 "P(S_t >= X, m_t >= B)", 
				 "P(M_t <= B)",
                 "P(M_t >= B)",
                 "P(m_t <= B)",
                 "P(m_t >= B)"), 
 		         S0=100, X, B, t=1, mu=0, sigma=1)
{
  Type = match.arg(Type)
  if (Type == "P(S_t <= X)") {
	# Bsp.: Shortfall-Risk für Plain-Vanilla-Call  Risiko, am Ende der Laufzeit aus dem Geld zu sein. Formal : P(S_T <= X)  
	pnorm((log(X/S0) - (mu - sigma^2/2)*t) / sqrt(sigma^2 * t), mean=0, sd=1);  
  } else if (Type == "P(S_t >= X)") {
	# Bsp.: Shortfall-Risk für Plain-Vanilla-Put  Risiko, am Ende der Laufzeit aus dem Geld zu sein. Formal : P(S_T >= X)  
	1 - (pnorm((log(X/S0) - (mu - sigma^2/2)*t) / sqrt(sigma^2 * t), mean=0, sd=1))  
  } else if (Type == "P(S_t >= X, m_t >= B)") {
    # Formel nach Poulsen (2004), S.7
    # P(S(t) >= X, m(t) >= B)
    pnorm((log(S0/X) + (mu - sigma^2/2)*t) / sqrt(t*sigma^2)) - (S0/B)^(1-2*mu/sigma^2) * pnorm((log(B^2/(S0*X)) + (mu - sigma^2/2)*t) / sqrt(t*sigma^2))
  } else if (Type == "P(M_t >= B)") {
    drift <- (mu - sigma^2/2)
    y <- log(B/S0)
    # Siehe Chuang (1996), P(M_t >= a), Gleichung (2.3)
	  p <- 1 - pnorm((y - drift * t) / (sigma*sqrt(t))) + exp(2 * drift * y/sigma^2) * pnorm((-y - drift * t) / (sigma*sqrt(t)))  
	  return(p)
  } else if (Type == "P(M_t <= B)") {
    drift <- (mu - sigma^2/2)
    y <- log(B/S0)
    # Siehe Chuang (1996), P(M_t >= a), Gleichung (2.3)
	  p <- 1 - pnorm((y - drift * t) / (sigma*sqrt(t))) + exp(2 * drift * y/sigma^2) * pnorm((-y - drift * t) / (sigma*sqrt(t)))  
	  return(1-p)
  }	  
  else if (Type == "P(m_t <= B)") {
	  # Verteilung des Minimums Min(t) bei Geometrischer Brownscher Bewegung:  P(Min(t) <= B)
    # Formel nach Ingersoll (1987), S.353 und Anpassung des Driftterms bei GBM  
    drift <- (mu - sigma^2/2)
	  y <- log(B/S0) # < 0
	  p <- pnorm((y - drift * t)/(sigma*sqrt(t))) + exp((2 * drift * y)/sigma^2) * pnorm((y + drift * t)/(sigma*sqrt(t)))
	  p  
  } else if (Type == "P(m_t >= B)") {
	  # P(m_t >= B) = 1 - P(m_t <= B)  
	  1 - calculateProbabilityGeometricBrownianMotion(Type == "P(m_t <= B)", S0, X, B, t, mu, sigma)  
  }
}

# alias
calcGBMProbability <- calculateProbabilityGeometricBrownianMotion

# P(M_t >= a), Gleichung (2.3)
cdf.M_t <- function(a, t=1, mu=0.05, sigma=1) {
  1 - pnorm((a - mu*t) / (sigma*sqrt(t))) + exp(2*mu*a/sigma^2) * pnorm((-a - mu*t) / (sigma*sqrt(t)))	
}

#
# P(m_t <= a), Gleichung (2.3), aber für das Minimum, a < 0
#
# 
cdf.m_t <- function(a, t=1, mu=0.05, sigma=1)
{
  1 - pnorm((-a + mu*t) / (sigma*sqrt(t))) + exp(2*mu*a/sigma^2) * pnorm((a + mu*t) / (sigma*sqrt(t)))	
}

# Bestimme gemeinsame Verteilungsfunktion für (M_t, B_t) für a>0
#
# P(M_t >= a, B_t <= z) vgl. Chuang (1996), Gleichung (2.1), S.82
#
joint.cdf.M_t.B_t <- function(a=0, z=0, t=1, mu=0.05, sigma=1)
{
  if (a < 0 ) stop("a muss groesser sein als 0")
  if (z <= a)
  {
  	exp(2*mu*a/sigma^2) * pnorm((z - 2*a - mu*t)/(sigma*sqrt(t)))
  }
  else # z > a
  {
    pnorm((z - mu*t)/(sigma*sqrt(t))) - pnorm((a - mu*t)/(sigma*sqrt(t))) + exp(2*mu*a/sigma^2) * pnorm((-a - mu*t)/(sigma*sqrt(t)))
  }
}

# Bestimme gemeinsame Verteilungsfunktion (m_t, B_t) (a < 0)
#
# P(m_t <= a, B_t >= z) vgl. Chuang (1996), Gleichung (2.1), S.82, für das Minimum
#
joint.cdf.m_t.B_t <- function(a=0, z=0, t=1, mu=0.05, sigma=1)
{
	if (a > 0) stop("a muss kleiner sein als 0")
	if (z >= a)
	{
	  exp(2*mu*a/sigma^2) * pnorm((-z + 2*a + mu*t)/(sigma*sqrt(t)))
	}
	else # z < a
	{
	  pnorm((-z + mu*t)/(sigma*sqrt(t))) - pnorm((-a + mu*t)/(sigma*sqrt(t))) + exp(2*mu*a/sigma^2) * pnorm((a + mu*t)/(sigma*sqrt(t)))
	}
}

# Bestimme P(m_t <= B) für die GBM
prob.barrier<-function(S0, B=60, t=2, mu=0.05, sigma=0.2)
{
  cdf.m_t(a=log(B/S0), t=t, mu=(mu-sigma^2/2), sigma=sigma)
}


# Gemeinsame Dichtefunktion P(S(t) >= x, Min(t) >= y) aus der Simulation bestimmen
joint.cdf.GBM.Min.Simulation <- function(x, y, S_t, Min_t)
{
   if (length(S_t) != length(Min_t)) stop("S_t und Min_t muessen die gleiche Laenge haben!")
   if (x < y) stop ("x>=y")
   
   mean(S_t >= x & Min_t >= y)
}

# 
# Formel nach Poulsen (2004), S.7
# P(S(t) >= k, Min(t) >= B)
joint.cdf.GBM.Min.Poulsen <- function(K, B, S0=100, T=1, mu=0, sigma=1)
{
  pnorm((log(S0/K) + (mu - sigma^2/2)*T) / sqrt(T*sigma^2)) - (S0/B)^(1-2*mu/sigma^2) * pnorm((log(B^2/(S0*K)) + (mu - sigma^2/2)*T) / sqrt(T*sigma^2))
}

# Verteilung des Minimums Min(t) bei Geometrischer Brownscher Bewegung:  P(Min(t) <= B)
# 
# Formel nach Ingersoll (1987), S.353 und Anpassung des Driftterms bei GBM
cdf.GBM.Min <- function(B, S0=100, T=1, mu=0, sigma=1)
{
	y <- log(B/S0) # < 0
	p <- pnorm((y - (mu - sigma^2/2) * T)/(sigma*sqrt(T))) + exp((2 * (mu - sigma^2/2) * y)/sigma^2) * pnorm((y + (mu - sigma^2/2) * T)/(sigma*sqrt(T)))
	p
}

# Formel (2.7) nach Chuang (1996), S.84 für das Maximum M_s und den Wert B_t (M_s,B_t): 
#
# P(M_s >= a, B_t <= z) mit s<t
#
joint.cdf.M_s.B_t1 <- function(a=110, z=100, s=0.6, t=1, mu=0.05, sigma=1)
{
	# Korrelation zwischen den Zeitpunkten s und t
  rho <- sqrt(s/t)
    
  # correlation matrix between s and t
  R   <- matrix(c(1  , rho, 
                 rho,  1), 2, 2)
  # 
  x <- (-a - mu*s)/(sigma*sqrt(s))
    
	p1 <- pmvnorm(upper=c(x, (z-2*a-mu*t)/(sigma*sqrt(t))), corr = R)
	p2 <- pmvnorm(upper=c((a - mu*s)/(sigma*sqrt(s)), (z    -mu*t)/(sigma*sqrt(t))), corr = R)
	
	exp(2*mu*a/sigma^2) * p1 + pnorm(z-mu*t/(sigma*sqrt(t))) - p2
}
