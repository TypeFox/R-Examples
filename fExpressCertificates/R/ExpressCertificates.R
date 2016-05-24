# Pricing of Express Certificate types
#
# Stefan Wilhelm

# Pricing of Easy-Express-Certificate
#
# Duplication:
# (1) Rückzahlungsbetrag S_0 (z.B. 120 EUR)
# (2) Cash-or-Nothing (short-put) mit Strike = B und Auszahlung von (S_0-B) : 
#     Ist das Underlying am Bewertungstag unter der Barriere, 
#     verliert man erstmal noch den Bonus.
# (3) Short-Put mit Strike = B
#
#
# @params S underlying price
# @params S0 Rückzahlungsbetrag bzw. Höchstbetrag (z.B. 120 EUR)
# @params B barrier (B <= S0)
# @params Time time to maturity in years
# @params r interest rate p.a. as 0.02 = 2%
# @params r_d continuous dividend yield p.a. as 0.02 = 2%
# @params sigma volatility p.a. as 0.18 = 18%
# @params ratio
EasyExpressCertificate<-function(S, S0, B, Time, r, r_d, sigma, ratio=1)
{
  # 1. Short Cash-or-Nothing-Put mit Strike X=B und Auszahlung (S0-B)
  con_short_put <- CashOrNothingOption("p", S, X=B, (S0-B), Time, r, b=r-r_d, sigma)
  price1 <- pmax(attr(con_short_put,"price"),0)
  
  # 2. Plain-Vanilla-Short-Put mit Strike B
  plain_short_put <- GBSOption(TypeFlag="p", S, X=B, Time, r, b=r-r_d, sigma)
  price2 <- pmax(attr(plain_short_put,"price"),0)
  
  (S0 * exp(-r*Time) - price1 - price2) * ratio
}

# payoff function g(S_T) at maturity
g <- function(S_T) { S_T }

#
# Bewertung Express-Classic-Zertifikat über Formel/numerische Integration  
#
# @param S
# @param X
# @param T
# @param K vector of cash rebates for early redemptions (length (n-1))
# @param g
ExpressCertificate.Classic <- function(S, X, T, K, g=function(S_T){S_T}, r, r_d, sigma, ratio=1)
{
	# Anzahl der Bewertungstage
	n <- length(T)
	
	# vorzeitige Stop-Wahrscheinlichkeiten p_0i
	p_0i <- calcRedemptionProbabilities(S=S, X=X, T=T, r=r, r_d=r_d, sigma=sigma)$stop.probs
    
	# Kovarianzmatrix
	Sigma <- sigma^2*outer(T,T,pmin)
	
	# Mittelwertvektor
	mu <- ((r - r_d) - sigma^2/2) * T
	
	# Unterer Trunkierungsvektor a<=x<=b der Renditen
	a <- rep(-Inf, n)
	
	# Oberer Trunkierungsvektor a<=x<=b
	b <- c(log(X) - log(S), Inf)
	
	# Randdichte der Renditen am Schluss
	f <- function(x)
	{
	  fx <- dtmvnorm.marginal(x, n=n, mean=mu, sigma=Sigma, lower=a, upper=b)
	  fx
	}
	
	# Definiere g(x)*f(x) fürs Integrieren E[g(X)] = int_{g(x)*f(x) dx} für die Variable i
	expectation <- function(x)
	{
		g(exp(x)*S) * f(x)
	}
	
	# Erwartungswert der Auszahlung am Laufzeitende E[g(S_T)]
	E_gST <- integrate(expectation, lower=-10, upper=10)$value
    
	# Preis = Erwartungswert aller diskontierten Zahlungen
	pi0 <- sum(c(K, E_gST) * p_0i * exp(-r * T))
	pi0
}

# Express-Zertifikat : Mehrere vorzeitige Bewertungstage, 
# am letzten Bewertungstag wie ein normaler Zero-Strike-Call
#
# @param S aktueller Stand des Basiswerts
# @param S0 Startkurs Basiswert (wenn Grenzen prozentual angegeben sind)
# @param B Barriere am Laufzeitende
# @param X Strike-Preise für n Bewertungstage (Vektor der Länge n)
# @param Cap Cap
# @param T Bewertungstage (Vektor der Länge n)
# @param Payoff Rückzahlungsvektor für vorzeitige Rückzahlungen (Vektor der Länge (n-1))
ExpressCertificate<-function(S, S0, X, Time, Payoff, r, r_d, sigma, ratio=1)
{
  # Anzahl der Bewertungstage
  n = length(Time)
  
  # Fürs Debugging: Vektor der diskontierten Rückzahlungspreise
  prices<-c()
  
  # Result-Vektor der Prices
  res<-c()
  
  # Express-Zertifikate am Laufzeitende: Time==0
  if (n==1 && Time==0)
  {
    zero_strike_call = GBSOption(TypeFlag="c", S=S, X=0, Time=0, r, b=r-r_d, sigma) 
    res=attr(zero_strike_call,"price")
  }
  # Express-Zertifikate vor Laufzeitende
  else
  {
    for (k in 1:length(S))
    {
      # unbedingte Wahrscheinlichkeit, dass das Zertifikat bis zum Bewertungstag j läuft : P(S_j>=X_j | S_i<=X_i für alle i < j)
      probs<-calculate.conditional.props(S=S[k], X=X, T=Time, r, r_d=r_d, sigma=sigma)$unconditional.probs
  
      for (i in 1:(n-1))
      {
        # Diskontierte Zahlung zum Zeitpunkt i gewichtet mit der Wahrscheinlichkeit des Stops zum Zeitpunkt i 
        prices[i]=probs[i]*(Payoff[i])/(1+r)^Time[i]
      }
  
      # letzter Bewertungstag : Bonuszertifikat gewichtet mit der Wahrscheinlichkeit des Stops zum Zeitpunkt n
      zero_strike_call = GBSOption(TypeFlag="c", S=S[k], X=0, Time=Time[n], r, b=r-r_d, sigma) 
      price_zero_strike_call = attr(zero_strike_call,"price")
      prices[n] <- probs[n] * price_zero_strike_call
      
      res[k] = sum(prices)
    }
  }
  res
}


# Express-Bonus-Zertifikat : Mehrere vorzeitige Bewertungstage, am letzten Bewertungstag wie ein normales Bonus-Zertifikat
#
# 1. Barriere B (=Pfadabhängigkeit greift nur am letzten Bewertungstag)
# 
# @param S aktueller Stand des Basiswerts
# @param S0 Startkurs Basiswert (wenn Grenzen prozentual angegeben sind)
# @param B Barriere am Laufzeitende
# @param X Strike-Preise für n Bewertungstage (Vektor der Länge n)
# @param Cap Cap
# @param T Bewertungstage (Vektor der Länge n)
# @param Payoff Rückzahlungsvektor für vorzeitige Rückzahlungen (Vektor der Länge (n-1))
ExpressBonusCertificate<-function(S, S0, B, X, Time, Payoff, r, r_d, sigma, ratio=1, DEBUG=FALSE)
{
  # Anzahl der Bewertungstage
  n = length(Time)
  
  # Fürs Debugging: Vektor der diskontierten Rückzahlungspreise
  prices<-c()
  
  # Result-Vektor der Prices
  res<-c()
  
  # Express-Zertifikate am Laufzeitende: Time==0
  if (n==1 && Time==0)
  {
    res=BonusCertificate(S=S, X=X[length(X)], B=B, Time=0, r, r_d, sigma, ratio=1)
  }
  # Express-Zertifikate vor Laufzeitende
  else
  {
    for (k in 1:length(S))
    {
      # unbedingte Wahrscheinlichkeit, dass das Zertifikat bis zum Bewertungstag j läuft : P(S_j>=X_j | S_i<=X_i für alle i < j)
      probs<-calculate.conditional.props(S=S[k], X=X, T=Time, r, r_d=r_d, sigma=sigma)$unconditional.probs
  
      for (i in 1:(n-1))
      {
        # Diskontierte Zahlung zum Zeitpunkt i gewichtet mit der Wahrscheinlichkeit des Stops zum Zeitpunkt i 
        prices[i]=probs[i]*(Payoff[i])/(1+r)^Time[i]
      }
  
      # letzter Bewertungstag : Bonuszertifikat gewichtet mit der Wahrscheinlichkeit des Stops zum Zeitpunkt n
      prices[n] <- probs[n]*BonusCertificate(S=S[k], X=X[n], B=B, Time=Time[n], r, r_d, sigma, ratio=1)
  
      #Fürs Debugging : 
      if (DEBUG)
      {
        res[k] = list(price=sum(prices), prices=prices, probs=probs)
      }
      else
      {
        res[k] = sum(prices)
      }
    }
  }
  res
}


# Express-Capped-Bonus-Zertifikat : Mehrere vorzeitige Bewertungstage, am letzten Bewertungstag wie ein normales Capped-Bonus-Zertifikat
#
# 1. Barriere B (=Pfadabhängigkeit greift nur am letzten Bewertungstag)
# 
# @param S aktueller Stand des Basiswerts
# @param S0 Startkurs Basiswert (wenn Grenzen prozentual angegeben sind)
# @param B Barriere am Laufzeitende
# @param X Strike-Preise für n Bewertungstage (Vektor der Länge n)
# @param Cap Cap
# @param T Bewertungstage (Vektor der Länge n)
# @param Payoff Rückzahlungsvektor für vorzeitige Rückzahlungen (Vektor der Länge (n-1))
ExpressCappedBonusCertificate<-function(S, S0, B, X, Cap, Time, Payoff, r, r_d, sigma, ratio=1, DEBUG=FALSE)
{
  # Anzahl der Bewertungstage
  n = length(Time)
  
  # Fürs Debugging: Vektor der diskontierten Rückzahlungspreise
  prices<-c()
  
  # Result-Vektor der Prices
  res<-c()
  
  # Express-Zertifikate am Laufzeitende: Time==0
  if (n==1 && Time==0)
  {
    res=CappedBonusCertificate(S=S, X=X[length(X)], B=B, Cap=Cap, Time=0, r, r_d, sigma, ratio)
  }
  # Express-Zertifikate vor Laufzeitende
  else
  {
    for (k in 1:length(S))
    {
      # unbedingte Wahrscheinlichkeit, dass das Zertifikat bis zum Bewertungstag j läuft : P(S_j>=X_j | S_i<=X_i für alle i < j)
      probs<-calculate.conditional.props(S=S[k], X=X, T=Time, r_f=r, r_d=r_d, sigma=sigma)$unconditional.probs
  
      for (i in 1:(n-1))
      {
        # Diskontierte Zahlung zum Zeitpunkt i gewichtet mit der Wahrscheinlichkeit des Stops zum Zeitpunkt i 
        prices[i]=probs[i]*(Payoff[i])/(1+r)^Time[i]
      }
  
      # letzter Bewertungstag : Bonuszertifikat gewichtet mit der Wahrscheinlichkeit des Stops zum Zeitpunkt n. Basierend auf heutiger Sicht, daher kein Abzinsen
      prices[n] <- probs[n] * CappedBonusCertificate(S[k], X=X[n], Cap=Cap, B=B, Time=Time[n], r=r, r_d=r_d, sigma=sigma, ratio=ratio)
      if (DEBUG) print(prices)
      res[k] = sum(prices)
    }
  }
  res
}

###############################################################################################################

# Hilfsfunktionen

###############################################################################################################

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
  #p = pnorm(log(X/S), mean=mu*T, sd = sqrt(sigma^2 * T))
  p = pnorm((log(X/S) - (mu - sigma^2/2)*T) / sqrt(sigma^2 * T), mean=0, sd=1);
  if (type == "c")
  {
    p
  }
  else
  {
   (1-p)
  }
}

# Hilfsfunktion zum Berechnen der Restlaufzeit in Jahren (RLZ) für eine Reihe von Daten
# Beispiel: RLZ(c("16.06.2008","16.06.2009","16.06.2010"))
#
# @param Vektor mit Datumsangabe
# @param Datumsformat, Standard ist "%d.%m.%Y"
RLZ<-function(dates, dateformat="%d.%m.%Y", start=NA)
{
  # Referenzdatum t0 ist aktuelles Datum
  t0 = ifelse(!is.na(start) & start != "", as.Date(start, format=dateformat), Sys.Date())
  
  # Restlaufzeiten in Jahren
  rlz = as.numeric(difftime(as.Date(dates ,format=dateformat), as.Date(t0)), units="days")/365
  rlz
}

annualizeYield<-function(yield, T)
{
  ifelse (T<=1, (yield/T), ((1+yield)^(1/T)-1))
}
