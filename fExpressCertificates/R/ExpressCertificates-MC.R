# Monte Carlo Methods to price derivative securities
#
#

# Monte-Carlo-Simulation für pfadunabhängige Classic-Expresszertifikate (Euler Scheme)
# 
# @param S aktueller Stand des Basiswerts
# @param X Strike-Preise für n Bewertungstage (Vektor der Länge (n-1))
# @param T Bewertungstage (Vektor der Länge n)
# @param K vector of cash rebates for early redemptions (length (n-1))/Rückzahlungsvektor für vorzeitige Rückzahlungen (Vektor der Länge (n-1))
# @param mc.steps Anzahl der Datenpunkte im Prozess
# @param mc.loops Anzahl der MC-Versuche
MonteCarlo.ExpressCertificate.Classic <- function(S, X, T, K, r, r_d, sigma, ratio=1, mc.steps=1000, mc.loops=20)
{
	# Anzahl der Bewertungstage
	n <- length(T)
	
	# M simulierte Preise, j=1..m
	prices <- numeric(mc.loops)
	
	# Preise am Laufzeitende
	S_T    <- rep(NA,mc.loops)
	
	# M Versuche, Zeitpunkt i=1..n, wann das Zertifikat stoppt bzw. wie lang es läuft.
	stops <- c()
	
	dt = max(T)/mc.steps
	
	for (j in 1:mc.loops)
	{
		# Simuliere den kompletten Prozess S(t)
		S_t <- GBM(S0=S, mu=r-r_d, sigma, T=max(T), N=mc.steps)
		early.redemption<-FALSE
		
		# Für alle Bewertungstage 1..(n-1)
		for (i in 1:(n-1))
		{
			# Aktienkurs am Bewertungstag T[i]--> S[T[i]/n+1]
			if (S_t[round(T[i]/dt+1)]>=X[i]) {
				prices[j] <- K[i]*exp(-r*T[i])
				stops[j] <- i
				early.redemption <- TRUE
				break
			}
		}  
		
		if (early.redemption==FALSE)
		{
			stops[j]  <- n
			S_T[j]    <- S_t[mc.steps]
			prices[j] <- S_t[mc.steps]*exp(-r*T[n])
		}
	}
	result=list(stops=stops, prices=prices, p=mean(prices), S_T=S_T)
	result
}

# conditional sampling from conditional densities f(x1|x0), f(x2|x1) ... f(xn|x_(n-1))
#
# @param S aktueller Stand des Basiswerts
# @param X Strike-Preise für n Bewertungstage (Vektor der Länge (n-1))
# @param T Bewertungstage (Vektor der Länge n)
# @param K vector of cash rebates for early redemptions (length (n-1))/Rückzahlungsvektor für vorzeitige Rückzahlungen (Vektor der Länge (n-1))
# @param conditional.random.generator can be "rnorm", "rnorm.halton", "rnorm.sobol"
Conditional.MonteCarlo.ExpressCertificate.Classic <- function(S, X, T, K, r, r_d, sigma, ratio=1, mc.loops=20, 
  conditional.random.generator="rnorm")
{
  # Anzahl der Bewertungstage
  n <- length(T)
	
  # M simulierte Preise, j=1..m
  prices <- numeric(mc.loops)
	
  # Preise zu den Bewertungstagen
  r_t    <- matrix(NA,mc.loops,n)
  S_t    <- matrix(NA,mc.loops,n)
	
  # M Versuche, Zeitpunkt i=1..n, wann das Zertifikat stoppt bzw. wie lang es läuft.
  stops<-c()
  
  mu  <- r-r_d
  
  ind <- rep(FALSE, mc.loops)
  	
  for (j in 1:n)
  {
    if (j == 1) 
    {
      dt <- T[1]
      x0 <- 0
    }
    else 
    {
      dt <- T[j] - T[j-1]
      x0 <- r_t[,j-1]
    }
    
	# Simuliere Renditen zum Zeitpunkt t_j
    if (conditional.random.generator == "rnorm")
    {
      r_t[,j] <- rnorm(n=mc.loops) * sqrt(sigma^2*dt) + x0 + (mu-sigma^2/2)*dt
    }
    else if (conditional.random.generator == "rnorm.halton")
    {
      # TODO : Ist das dann noch low-discrepancy in 2 Dimensionen ?
      r_t[,j] <- sample(rnorm.halton(n=mc.loops, dimension=1, init=TRUE)) * sqrt(sigma^2*dt) + x0 + (mu-sigma^2/2)*dt
    }
    else if (conditional.random.generator == "rnorm.sobol")
    {
      r_t[,j] <- rnorm.sobol(n=mc.loops, dimension=1, init=TRUE, scrambling=TRUE) * sqrt(sigma^2*dt) + x0 + (mu-sigma^2/2)*dt
    }
    else
    {
      r_t[,j] <- rnorm(n=mc.loops, mean=x0+(mu-sigma^2/2)*dt, sd=sqrt(sigma^2*dt))
    }
    
    S_t[,j] <- S * exp(r_t[,j])
        
    if (j < n)
    {
      ind2         = !ind & S_t[,j] > X[j]
    }
    else
    {
      ind2         = !ind
    }
    stops[ind2] = j 
    #table(stops)
    if (j < n)
    {
      prices[ind2]  =   K[j] * exp(-r*T[j]) 
    }
    else
    {
      prices[ind2]  = S_t[ind2,j] * exp(-r*T[j]) 
    }
    ind = ind | ind2
  }
  result=list(stops=stops, prices=prices, p=mean(prices), S_T=S_t[,n])
  result
}

# (Naive/Euler Scheme) Monte-Carlo-Bewertung für pfadunabhängige Express-Zertifikate 
# mit variabler Auszahlungsstruktur payoff.function am Laufzeitende
#
# 
# @param S aktueller Stand des Basiswerts
# @param X Strike-Preise für n Bewertungstage (Vektor der Länge n)
# @param T Bewertungstage (Vektor der Länge n)
# @param K Rückzahlungsvektor für vorzeitige Rückzahlungen (Vektor der Länge (n-1))
# @param B Barriere am Laufzeitende
# @param mc.steps Anzahl der Datenpunkte im Prozess
# @param mc.loops Anzahl der MC-Versuche
MonteCarlo.ExpressCertificate <- function(S, X, T, K, B, r, r_d, sigma, mc.steps=1000, mc.loops=20, 
  payoff.function)
{
  # Anzahl der Bewertungstage
  n <- length(T)
  
  # M simulierte Preise, j=1..m
  prices <- c()
  
  # M Versuche, Zeitpunkt i=1..n, wann das Zertifikat stoppt bzw. wie lang es läuft.
  stops<-c()
  
  for (j in 1:mc.loops)
  {
    # Simuliere den kompletten Prozess S(t)
    dt  <- max(T)/mc.steps
    S_t <- GBM(S0=S, mu=r-r_d, sigma, T=max(T), N=mc.steps)
    early.redemption<-FALSE
    
    # Für alle vorzeitigen Bewertungstage 1..(n-1)
    for (i in 1:(n-1))
    {
      # Aktienkurs am Bewertungstag T[i]--> S[T[i]/n+1]
      if (S_t[round(T[i]/dt+1)] >= X[i])
      {
        prices[j] <- exp(-r * T[i]) * K[i] 
        stops[j]  <- i
        early.redemption=TRUE
        break
      }
    }  
    
    # Payoff am Laufzeitende abhängig vom realisierten Pfad S_t ausrechnen
    if (early.redemption==FALSE)
    {
      stops[j]   <- n
      prices[j]  <- exp(-r * T[n]) * payoff.function(S_t)
    }
  }
  result=list(stops=stops, prices=prices, price=mean(prices), prob.stops=table(stops)/mc.loops)
  result
}

# Returns redemption index i, such that S(t_i) > X(t_i) and all S(t_j) < X(t_j) (j < i).
#
# @param S Vektor mit n Bewertungstagen (evtl. mehr Spalten)
# @param n Anzahl der Bewertungstage
# @param X Vektor der Strikepreise (Länge (n-1))
#
# @return Zeitpunkt der Rückzahlung 1..n
getRedemptionTime = function(S, n, X) {
	earlyredemption <- FALSE
	redemptionTime <- n
	for (i in 1:(n-1))
	{
		# Aktienkurs am Bewertungstag T[i]--> S[T[i]/n+1]
		earlyredemption <- S[i] >= X[i]
		if (earlyredemption) {
			redemptionTime <- i
			break
		}
	}
	redemptionTime
}

getRedemptionTime2 = function(S, n, X) 
{
	ind = 1:(n-1)
	if (any(S[ind] >= X[ind])) {
		return(min(which(S[ind] >= X[ind])))
	} else {
		return(n)
	}
}

# determine redemption times for a matrix of underlying prices S (N x n).
# should be much faster than apply(S, 1, getRedemptionTime).
#
# @param S (N x n) matrix
# @param n
# @param X (n x 1) vector (call level for t_i is then X[i]) or (N x n) matrix (call levels for t_i and case j=1...N are then X[j,i]) 
# @return (N x 1) mit Werten aus 1..n
getRedemptionTimesForMatrix <- function (S, n, X) {
   N <- nrow(S)
   redemptionTime <- rep(n, N)       # redemption at t_n by default
   for (i in 1:(n - 1)) {
     if (length(dim(X)) == 2) {
       ind_earlyredemption <- S[,i] >= X[,i]
     } else {
       ind_earlyredemption <- S[,i] >= X[i]
     }
     redemptionTime[ind_earlyredemption] <- pmin(redemptionTime[ind_earlyredemption], i)
   }
   return(redemptionTime)
}

# Bestimmt die Auszahlung für einen konkreten Rückzahlungszeitpunkt i = 1..n
# @param i Rückzahlungszeitpunkt i = 1..n
# @param S Vektor der Länge n mit den Kursen an den n Bewertungstagen
# @param m_t Minimum bis zum Zeitpunt T
payoffExpressClassic <- function(i, n, S, m, K) {
 if (i < n) 
 	return(K[i])
 else 
	return(S[n])
}


# Payoff für ein Express-Bonus: Typ 1
# Wenn min(S_T) <= B, dann Auszahlung S_T/S0*100, 
# sonst K_T
# 
# Bestimmt die Auszahlung für einen konkreten Rückzahlungszeitpunkt i = 1..n
# @param i Rückzahlungszeitpunkt i = 1..n
# @param S Vektor der Länge n mit den Kursen an den n Bewertungstagen
# @param m_t Minimum bis zum Zeitpunt T
payoffExpressCappedBonusType1 <- function(i, n, S, m, K, B) {
 if (i < n) {
   return(K[i])
 } else {
   if (m[n] <= B) {
      return(S[n])
   } else {
      K[n]
   }
 }
}


# Payoff für ein Express-Bonus: Typ 2
# Wenn min(S_T) <= B, dann Auszahlung S_T/S0*100, 
# sonst max(K_T, S_T/S0*100)
# 
# Bestimmt die Auszahlung für einen konkreten Rückzahlungszeitpunkt i = 1..n
# @param i Rückzahlungszeitpunkt i = 1..n
# @param S Vektor der Länge n mit den Kursen an den n Bewertungstagen
# @param m_t Minimum bis zum Zeitpunt T
payoffExpressBonusType1 <- function(i, n, S, m, K, B) {
 if (i < n) {
   return(K[i])
 } else {
   if (m[n] <= B) {
      return(S[n])
   } else {
      max(K[n],S[n])
   }
 }
}

# Bestimmt die Auszahlung für einen konkreten Rückzahlungszeitpunkt i = 1..n
# @param i Rückzahlungszeitpunkt i = 1..n
# @param S Vektor der Länge n mit den Kursen S_t an den n Bewertungstagen
# @param m Vektor der Länge n mit den Minima m_t an den n Bewertungstagen
# TO BE DONE Kuponereignisse falls nicht vorzeitig zurückgezahlt, aber Barriere noch intakt
payoffExpressML0AN5 <- function(i, n, S, m, K, B, S0) {
 # Auszahlungsvektor inklusive Kuponzahlungen
 kupon      <- rep(0, n)
 redemption <- rep(0, n)
 
 # Schauen, ob wir Kuponzahlungen gekriegt hätten
 j = 1
 while (j < i) {
   if (m[j] <= B) {
     kupon[j] = 0
   }
   else {
     kupon[j] = 1 * S[j] / 1000
   }
   j  = j + 1
 }
 
 if (i < n) {
   # vorzeitige Kündigung
   redemption[i]     = (100 + 1 * S[i] / 1000)
 } else {
   # Am letzten Bewertungstag
   if (m[n] <= B) {
       redemption[n] = 100 * S[n] / S0       # prozentuale Rückzahlung bei Barrierenereignis
   } else {
      redemption[n]  = 100 + 1 * S[n] / 1000 # Rückzahlung plus Kupon wenn Barriere nicht getroffen wurde
   }
 }
 return (kupon + redemption)
}

# vgl. bessere Methode simExpressPriceMVN in express-simulation.R
# Probleme: payoff function cannot deal with memory mechanism
#
# 1. simulate joint vector of prices and minimum (S(t_1),...,S(t_n); m(t_1),...,m(t_n))' at valuation dates using Euler scheme
# 2. 
# @param payoffFunction is a function  payoffFunction(i, n=n, S=as.numeric(S[j,]), m=S[j,minimumFids], K,...)
#                       that determines the payoff at time t_i
SimulateGenericExpressCertificate <- function(S, X, K, T, r, r_d, sigma, 
  mc.loops=10000, mc.steps=1000, payoffFunction=payoffExpressClassic, ...)
{
	# 0. Sortiere Zeiten in der Vergangenheit aus
	# Schauen, ob Bewertungstage bereits vorbei sind
	ind <- which(T<0)
	# Wenn Bewertungstag bereits vorbei, dann ausschliessen
	if (length(ind)>0)
	{
		T      <- T[-ind]
		X      <- X[-ind]
		K      <- K[-ind]
	}
	
	# number of valuation dates
	n <- length(T)
	
	# 1. simulate joint vector of prices and minimum (S(t_1),...,S(t_n); m(t_1),...,m(t_n))' 
  #    at valuation dates using Euler scheme
  # see also: simPricesAndMinimumGBM
	E <- simPricesAndMinimumFromGBM2(N=mc.loops, S=S, T=T, mu=r-r_d, sigma=sigma, mc.steps)
	
	# 2. determine vector of redemption times i \in 1..n
	redemptionTimes  <- apply(X=E, MARGIN=1, FUN=getRedemptionTime, n = n, X)
    
  # column names for minimum
  minimumFids <- paste("m",1:n, sep="")
	
	# 3. discounted payoffs for redemption in t_i depending on prices and minimum
  #    (S(t_1),...,S(t_n); m(t_1),...,m(t_n))'
	prices <- sapply(seq(along=redemptionTimes), function(j, n, S, T, K, r, ...) {
        i = redemptionTimes[j] 
        exp(-r * T[i]) * payoffFunction(i, n=n, S=as.numeric(S[j,]), m=S[j,minimumFids], K, ...)
                        }, n, S=E, T, K, r, ...)
	
	l <- list(
      prices          = prices, 
			price           = mean(prices), 
			redemptionTimes = redemptionTimes, 
			n               = n,
			S               = S,
			X               = X,
			K               = K,
			T               = T)
	class(l) = "express.certificate"
	return(l)
}

SimulateExpressClassicCertificate <- function(S, X, K, T, r, r_d, sigma, mc.loops=10000, mc.steps=1000)
{
	SimulateGenericExpressCertificate(S, X, K, T, r, r_d, sigma, mc.loops, mc.steps, 
    payoffFunction=payoffExpressClassic)
}

# Payoff ist
# Wenn min(S_t) <= B, dann Auszahlung S_t/S0*100, 
# sonst max(K_n, S_t/S0*100)
SimulateExpressBonusCertificate <- function(S, X, B, K, T, r, r_d, sigma, mc.loops=10000, mc.steps=1000, barrierHit=FALSE)
{
  if (barrierHit) {
    SimulateGenericExpressCertificate(S, X, K, T, r, r_d, sigma, mc.loops, mc.steps, payoffFunction=payoffExpressClassic)
  } else {
    SimulateGenericExpressCertificate(S, X, K, T, r, r_d, sigma, mc.loops, mc.steps, payoffFunction=payoffExpressBonusType1, B=B)
  }
}

SimulateExpressCertificateML0AN5 <- function(S, X, K, T, r, r_d, sigma, mc.loops=10000, mc.steps=1000, B, S0, barrierHit=FALSE)
{
	# # 0. Sortiere Zeiten in der Vergangenheit aus
	# Schauen, ob Bewertungstage bereits vorbei sind
	ind <- which(T<0)
	# Wenn Bewertungstag bereits vorbei, dann ausschliessen
	if (length(ind)>0)
	{
		T      <- T[-ind]
		X      <- X[-ind]
		K      <- K[-ind]
	}
	n <- length(T)
	
	# 1. simulate joint vector of prices and minimum (S(t_1),...,S(t_n); m(t_1),...,m(t_n))' 
  #    at valuation dates using Euler scheme
  # see also: simPricesAndMinimumGBM
	E <- simPricesAndMinimumFromGBM2(N=mc.loops, S=S, T=T, mu=r-r_d, sigma=sigma,  mc.steps)
	
	# 2. Redemptionzeitpunkte 1..n bestimmen
	redemptionTimes <- apply(X=E, MARGIN=1, FUN=getRedemptionTime, n = n, X)
    
  # Spaltenbezeichnungen fürs Minimum und die Kurse
  minimumFids <- paste("m",1:n, sep="")
  kursFids    <- paste("S",1:n, sep="")
    
  # 3. Diskontierte Preise für ein Classic-Express bestimmen in Abhängigkeit vom Redemption-Zeitpunkt i und Kursverlauf, Verlauf des Minimums m_t
	prices <- sapply(seq(along=redemptionTimes), function(j, n, S, T, K, r, B, S0) {
     i = redemptionTimes[j] 
     sum(exp(-r * T) * payoffExpressML0AN5(i, n=n, S=as.numeric(S[j,kursFids]), m=as.numeric(S[j,minimumFids]), K,  B, S0))
     }, n, S=E, T, K, r, B, S0)
	
	l <- list(
      prices          = prices, 
			price           = mean(prices), 
			redemptionTimes = redemptionTimes, 
			n               = n,
			S               = S,
			X               = X,
			K               = K,
			T               = T,
      B               = B,
      barrierHit      = barrierHit, 
      probability.barrier.hit.GBM      = mean(E[,minimumFids[n]]<=B, na.rm=TRUE),                     # "normale" Barrierenwahrscheinlichkeit einer GBM
      probability.barrier.hit.maturity = mean(E[redemptionTimes == n,minimumFids[n]]<=B, na.rm=TRUE)) # Barrierenwahrscheinlichkeit nur für die Pfade, die bis zum Ende überleben
            
	class(l) <- "express.certificate"
	return(l)
}

# print summary information for class "express.certificate"
print.express.certificate <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
	# TO BE DONE : Berechnung Mittelwert der MC-Preise/Standardabweichung/prozentuale Stop-Wahrscheinlichkeiten usw.
	#prob.stops            = table(stops)/M
	#prob.early.redemption = 1-(table(stops)/M)[n]
	#barrierHit            = barrierHit
	#prob.barrier.hit=mean(barrierHit, na.rm=TRUE)
	
	cat("\nExpress Certificate:\n\n", sep = "")
	cat("Price:\n")
	print.default(format(x$price, digits = digits), print.gap = 2, quote = FALSE)
	cat("\n")
	
	cat("Current underlying price:\n")
	print.default(format(x$S, digits = digits), print.gap = 2, quote = FALSE)
	cat("\n")
	
	# n = Anzahl der Bewertungstage
	n <- x$n
	
	# Stammdaten tabellarisch anzeigen
	t1 <- rbind(T = x$T, 
			X = c(x$X[-n], NA),
			K = c(x$K[-n], NA))
	colnames(t1)      = paste("Valuation date",1:n)
	dimnames(t1)[[1]] = c("Years T", "Strike price X", "Cash rebate K")
	print(t1)
	cat("\n")
	
	# prozentuale Rückzahlungswahrscheinlichkeiten
	t <- table(factor(x$redemptionTimes, levels=1:n))/length(x$redemptionTimes)
	cat("Redemption times probabilities:\n")
	print(t)
	cat("\n")
	
	cat("Probability of early redemption:\n")
	probability.early.redemption = 1-t[n]
	print.default(format(probability.early.redemption, digits = digits), print.gap = 2, quote = FALSE)
    
  if (!is.null(x$probability.barrier.hit.GBM) && !is.null(x$B)) {
   cat("Barrier:\n")
   print.default(x$B, print.gap = 2, quote = FALSE)
   cat("\n")
      
   cat("Probability of barrier event until maturity (all paths):\n")
   print.default(format(x$probability.barrier.hit.GBM, digits = digits), print.gap = 2, quote = FALSE)
   cat("\n")
      
   cat("Probability of barrier event conditioned on final redemption:\n")
   print.default(format(x$probability.barrier.hit.maturity, digits = digits), print.gap = 2, quote = FALSE)
   cat("\n")
	}
    
	cat("\n")
	invisible(x)
}


# examples/usage
if (FALSE) {
library(fExpressCertificates)
  
# 1. simulate joint vector of prices and minimum (S(t_1),...,S(t_n); m(t_1),...,m(t_n))' at valuation dates
# (a) using Euler scheme	
E <- simPricesAndMinimumFromGBM2(N=10000, S=100, T=c(0.5, 0.7, 1), mu=0.05, sigma=0.4, mc.steps=1000)
str(E)  # 10000 x 6

# (b) using MVN and minimum of Brownian Bridge	
E2 <- simPricesAndMinimumFromGBM(N=10000, S=100, T=c(0.5, 0.7, 1), mu=0.05, sigma=0.4)
str(E2) # 10000 x 4

# double check if E and E2 produce similar densities
plot(density(E[,"m3"], to=100), col="blue")
lines(density(E2[,"m3"], to=100), col="red")
	
# 2. Redemptionzeitpunkte 1..n bestimmen
X <- c(100, 100)
redemptionTimes  <- apply(X=E, MARGIN=1, FUN=getRedemptionTime, n = 3, X)
redemptionTimes2 <- getRedemptionTimesForMatrix(S=E, n=3, X)
all.equal(redemptionTimes, redemptionTimes2)
	
# Diskontierte Preise für ein Classic-Express bestimmen in Abhängigkeit vom Redemption-Zeitpunkt i
T <- c(0.5, 0.7, 1)
K <- c(120,130)
n <- 3
r <- 0.05

# get discounted payoffs (N x 1)	
prices <- sapply(seq(along=redemptionTimes), function(j) {
	i = redemptionTimes[j]
	if (i < n) 
		return(exp(-r * T[i]) * K[i])
	else 
		return(exp(-r * T[i]) * E[j,n])
})
mean(prices)

	
# künstliches Beispiel profilen
Rprof("sim1.out")
p <- SimulateExpressClassicCertificate(S=100, X=c(100,100), K=c(120,130), T=c(0.5, 0.7, 1), r=0.05, r_d=0, sigma=0.4, mc.loops=10000, mc.steps=1000)
Rprof(NULL)
summaryRprof("sim1.out")
table(p$redemptionTimes)
p$price
	
# Beispiel CB7AXR auf DTEGn.DE am 10.12.2009
p <- SimulateExpressBonusCertificate(S=10.4/12.10*100, X=c(100,100,100), B=7/12.1*100, K=c(134, 142.5, 151), 
			T=.RLZ(c("16.12.2009","17.06.2010","17.12.2010")), r=0.01, r_d=0, sigma=0.23, mc.loops=10000, mc.steps=1000)
p
    
# andere Methode
p2 <- SimulateExpressBonusCertificate2(S=10.4/12.10*100, X=c(100,100,100), B=7/12.1*100, K=c(134, 142.5, 151), 
			T=.RLZ(c("16.12.2009","17.06.2010","17.12.2010")), r=0.01, r_d=0, sigma=0.23, mc.loops=10000, mc.steps=1000)
p2
    
plot(density(p$prices))
lines(density(p2$prices), col="red")
	
# Beispiel BLB248 auf .STOXX50E am 10.12.2009
# Barriere wurde schon getroffen, daher wie Express Classic
p <- SimulateExpressClassicCertificate(S=2873.08 / 4401.68 * 100, X = c(100,100), K=c(130,140), T=.RLZ(c("06.09.2010","06.09.2011")), 
			r=0.01, r_d=0.02, sigma=0.3, mc.loops=10000, mc.steps=1000)
    table(p$redemptionTimes)
    p
plot(density(p$prices))

################################################################################

# Beispiel "DAX Kupon Zertifikat" ML0AN5 auf DAX zur Emission am 30.04.2008
# Variable Kuponzahlungen

# Startkurs zur Emission
S0 <- 6948.82

# @param S (N x n)
# @param m (N x n) Matrix der running Minima (m_1,m_2,m_3,m_4)
payoffML0AN5 <- function(S, m, X, K, B=50, S0=6948.82) {
  # Anzahl der Bewertungstage
  n <- ncol(S)
  N <- nrow(S)
  # lazy evaluation for B and S0
  S0 <- 6948.82
  B <- 50
  
  redemptionIndex <- getRedemptionTimesForMatrix(S, n, X)
    
  # Redemption Payoff Matrix (N x n)
  payoffMatrix <- matrix(0, N, n)
  for (i in 1:(n-1)) {
    # Bestimme variable Kuponzahlungen = Index_t / 1000 zum Zeitpunkt t_i
    coupon <- ifelse(m[,i] > B, (S[,i] * S0 / 100) / 1000, 0)
    
    # early redemption is always with coupon payments.
    ind_earlyredemption <- S[,i] >= X[i] & redemptionIndex == i
    payoffMatrix[ind_earlyredemption,i] <- 100 + (S[ind_earlyredemption,i] * S0 / 100) / 1000
    
    ind_coupon <- S[,i] < X[i] & redemptionIndex > i
    payoffMatrix[ind_coupon,i] <- coupon[ind_coupon]
  }
  ind_finalredemption <- redemptionIndex == n
  # Bestimme variable Kuponzahlungen Zeitpunkt t_n
  coupon <- ifelse(m[ind_finalredemption,n] > B, S[,n] * S0 / 100 / 1000, 0)
  payoffMatrix[ind_finalredemption, n] <- ifelse(m[ind_finalredemption, n] > B, 
    100 + coupon, S[ind_finalredemption, n])
    
  return(list(redemptionIndex=redemptionIndex, payoffMatrix=payoffMatrix))
}
	

# heutiger Kurs
S <- 6948.82 / S0*100

# 
T <- .RLZ(c("08.05.2009","07.05.2010","09.05.2011","08.05.2012","07.05.2013"), start="30.04.2008")

# Call Level 120% vom Startniveau
X <- c(120, 120, 120, 120) * S0 / 100
    
# Keine festen Auszahlungen, die Höhe der vorzeitigen Auszahlungen hängt von S_t ab!
K <- c(0, 0, 0)

# Barriere 50%
B <- 50 * S0 / 100

p <- SimulateExpressCertificateML0AN5(S=6948.82, X=X, K=K, T=T, r=0.01, r_d=0, sigma=0.24, B=B, S0=6948.82)
p

p2 <- simExpressPriceMVN(S=100, m=Inf, X=c(120, 120, 120, 120), K, B=50, T=T, r=0.01, r_d=0, sigma=0.24, 
  mc.loops = 100000, payoffFunction = payoffML0AN5)
    
# TODO: Formel-Notation, weil z.B. bei ML0AN5 der Payoff an den vorzeitigen Ausübungstagen hängt vom Indexstand ab!!!
}



 

