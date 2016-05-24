# Bestimmt die bedingten Wahrscheinlichkeiten P(S[i]>=X[i] | S[i-1]<=X[i-1]) 
# über die Brownsche Bewegung als Stop-Wahrscheinlichkeiten eines Express-Zertifikats
#
# @param S aktueller Aktienkurs
# @param X Vektor der Strike-Preise (Länge n)
# @param T Vektor der Restlaufzeiten (Länge n)
# @param r risikoloser Zinssatz
# @param r_d Dividendenrendite
# @param sigma Volatilität in % p.a.
#
# @return 
calculate.conditional.props <- function(S, X, T, r_f, r_d, sigma)
{
  mu <- r_f - r_d
  
  # Log-Renditen ausgehend vom heutigen Kurs
  r  <- log(X/S)
  
  # Erwartungswert der Log-Renditen für alle vorgegebenen Zeitpunkte T
  e_r <- (mu - sigma^2/2)*T
  
  # conditional Stop-Wahrscheinlichkeit vom Zeitpunkt i zum Zeitpunkt j: p_ij = P(S[j]>=X[j] | S[i]>=X[i])
  probs <- c()
  
  # Stop-Wahrscheinlichkeiten zum Zeitpunkt j : p_0j = P(S[j]>=X[j] & S[i=1..(j-1)]<=X[i=1..(j-1)])
  probs2 <- c()
  
  # Wahrscheinlichkeit P(R[1] >= r[1])
  probs[1]  <- 1-pnorm(r[1], mean=e_r[1], sd=sigma * sqrt(T[1]))
  probs2[1] <- probs[1]
  
  if (length(X) > 1)
  {
    for (i in 2:length(X))
    {
      # Bestimme die gemeinsame Wahrscheinlichkeit P(R_{i-1}<=r[i-1] und R_i >= r[i]) aus der Verteilungsfunktion der multivariaten NV.
    
      # Gemeinsame Normalverteilung der Log-Renditen r[i] und r[i-1].
      # theoretische Korrelation zwischen 2 Zeitpunkten s und t: rho(s,t) = min(s,t)*sigma^2/sqrt(s*sigma^2*t*sigma^2) = min(s,t) / sqrt(s*t)
      # Kovarianzmatrix zwischen beiden Zeitpunkten
      Sigma=matrix(c(  sigma^2 * T[i-1]                      ,  min(T[i-1],T[i]) * sigma^2,
                       min(T[i-1],T[i])  * sigma^2           ,  sigma^2 * T[i])
                     , 2, 2)
                     
      # Erwartungswert der Renditen                 
      mean=c(e_r[i-1],e_r[i])
                
      # Bestimme Bedingte Stop-Wahrscheinlichkeit p_{ij} : P(R_i >= r[i] | R[i-1] <= r[i-1]) = P(R[i-1] <= r[i-1] & R_i >= r[i])
      probs[i]  = pmvnorm(lower=c(-Inf,r[i]), upper=c(r[i-1],Inf), mean=mean, sigma=Sigma) / pmvnorm(lower=c(-Inf,-Inf), upper=c(r[i-1],Inf), mean=mean, sigma=Sigma)
    
      # unbedingten Wahrscheinlichkeiten p_{0j}
      probs2[i] = prod(1-probs[1:(i-1)]) * probs[i]
    }
  }
  
  probs2[length(T)] = 1-sum(probs2[1:(length(T)-1)])
  
  res=list(conditional.probs=probs,unconditional.probs=probs2)
  res
}


# Bestimmt die Stop-Wahrscheinlichkeiten P(S[i]>=X[i], S[1]<=X[1],...,S[i-1]<=X[i-1]) 
# über die Geometrische Brownsche Bewegung
#
# @param S aktueller Aktienkurs
# @param X Vektor der Strike-Preise (Länge n)
# @param T Vektor der Restlaufzeiten (Länge n)
# @param r risikoloser Zinssatz
# @param r_d Dividendenrendite
# @param sigma Volatilität in % p.a.
calcRedemptionProbabilities <- function(S, X, T, r, r_d, sigma) {
  mu <- r - r_d
  
  # Log-Renditen ausgehend vom heutigen Kurs
  r <- log(X/S)
  
  # Erwartungswert der Log-Renditen für alle vorgegebenen Zeitpunkte T
  e_r <- (mu - sigma^2/2)*T
  
  # Kovarianzmatrix
  # Gemeinsame Normalverteilung der Log-Renditen r[1] ... r[n]. Jeweils 2 Zeitpunkte r[i] und r[i-1] haben die Kovarianz = min(s,t) * sigma^2
  Sigma <- outer(T, T, pmin) * sigma^2
   
  # Stop-Wahrscheinlichkeiten zum Zeitpunkt j : p_0j = P(S[j]>=X[j] , S[i=1..(j-1)]<=X[i=1..(j-1)])
  probs <- c()
  
  # Wahrscheinlichkeit P(R[1] >= r[1])
  probs[1] <- 1-pnorm(r[1], mean=e_r[1], sd=sigma * sqrt(T[1]))
  
  for (i in 2:length(T)) {
    # P(S[2]>X[2],S[1]<X[1])
    if (i>1) {
      upper <- c(r[1:(i-1)],Inf)
    }
    else
    {
      upper = Inf
    }
    # Am letzten Termin
    if (i == length(T)) {
      lower <- rep(-Inf, i)
    }
    else {
      lower <- c(rep(-Inf, (i-1)), r[i])
    }
    probs[i] <- pmvnorm(lower=lower, upper=upper, mean=e_r[1:i], sigma=Sigma[1:i,1:i])
  }
   
  res <- list(stop.probs=probs)
  res
}
calcStopProbabilities <- calcRedemptionProbabilities

# Bestimmt die Stop-Wahrscheinlichkeiten P(S[i]>=X[i], S[1]<=X[1],...,S[i-1]<=X[i-1]) 
# über Monte Carlo Simulation
simRedemptionProbabilities <- function(S, X, T, r, r_d, sigma, mc.steps=1000, mc.loops=20)
{
	# Anzahl der Bewertungstage
	n <- length(T)
	
	# M Versuche, Zeitpunkt i=1..n, wann das Zertifikat stoppt bzw. wie lang es läuft.
	stops <- c()
	
	# dt
	dt <- max(T)/mc.steps
	
	# Index values
	ind <- round(T/dt+1)
	
	for (j in 1:mc.loops) {
		# Simuliere den kompletten Prozess S(t)
		S_t <- GBM(S0=S, mu=r-r_d, sigma, T=max(T), N=mc.steps)
		early.redemption <- FALSE
		
		# Für alle Bewertungstage 1..(n-1)
		for (i in 1:(n-1)) {
			# Aktienkurs am Bewertungstag T[i]--> S[T[i]/n+1]
			if (S_t[ind[i]] >= X[i])
			{
				stops[j] <-i
				early.redemption <- TRUE
				break
			}
		}  
		if (early.redemption==FALSE) {
		  stops[j] <- n
		}
	}
	result <- list(stops=stops)
	result
}
MonteCarloStopProbabilities <- simRedemptionProbabilities

