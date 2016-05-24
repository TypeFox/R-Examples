# Simulate simulate joint vector of prices and minimum (S(t_1),...,S(t_n); m(t_1),...,m(t_n))' 
# at valuation dates using various methods, e.g. 
# (1) Euler scheme
# (2) Multivariate normal distribution of returns at (t_1,...,t_n) and minimum of Brownian bridges
#     between t_{i-1} and t_i.
# (3)

# Simuliert mittels Monte Carlo (Euler Scheme) eine Matrix von Kursen aus der
# (unrestricted!!) Geometrischen Brownschen Bewegung an den Bewertungstagen eines Expresszertifikats 
# plus das Minimum mit folgender Struktur
#
# X_i = (S(t_1),...,S(t_n); m(t_1),...,m(t_n))'
#
# Diese Matrix kann dann als Grundlage für jedwede Auszahlungsstruktur (Express-Classic, Express-Bonus usw.) hergenommen werden.
# Das Minimum des Prozesses wird aus den Zwischenschritten (mc.steps) ermittelt. Wegen der vielen Zwischenschritte kann die Simulation sehr lange dauern.
#
# @param T valuation dates (length n)
# @param mu
# @param sigma
# @param mc.loops
# returns a matrix (mc.loops x 2*n) with (S(t_1),...,S(t_n); m(t_1),...,m(t_n))'
simPricesAndMinimumFromGBM2 <- function(N=10000, S=100, T, mu, sigma, 
  mc.steps=1000)
{
	n <- length(T)
	
	# Zeitintervall
	dt <- max(T)/mc.steps
	
	# vector of grid times
	t <- seq(0, T[n] , by=dt)
	
	# start value
	S0 <- S
	
	# matrix S with price paths
	S <- matrix(NA, N, mc.steps + 1)
	
	# prices at n valuation dates n (t_1,...,t_n): (S(t_1),...,S(t_n))'
	S_t <- matrix(NA, N, n)
	colnames(S_t) <- paste("S",1:n, sep='')
    
  # running minimum at n valuation dates i = 1...n
	m_t <- matrix(NA, N, n)
	colnames(m_t) <- paste("m",1:n, sep='')
	
	# Kurse simulieren
	for (j in 1:N)
	{
		S[j,] <- GBM(S0 = S0, mu=mu, sigma=sigma, T=T[n], N=mc.steps)
	}
	
	# Spalten-Indizes der Verfallstermine im Array S
	it  <- round(T/dt+1)
	
	# kumulierte Minima bestimmen Matrix(mc.loops x (mc.steps + 1))
	m <- t(apply(S, MARGIN=1, FUN=cummin))
	
	# Verteilung der Kurse S(t_i) zu den n Verfallsterminen, wenn der Pfad bis zum Ende durchläuft 
	for (i in 1:n)
	{
	  S_t[,i] <- S[,it[i]]
	  
	  # Minimum pro Pfad bis it[i] bestimmen
    m_t[,i] <- m[,it[i]]
	}
    
	# Rückgabe der Matrix (S_{t1}, S_{t2}, ..., S_{tn}, m_{tn})
	return(data.frame(cbind(S_t, m_t)))
}

# Simuliere Samples aus der (n+1)-dimensionalen gemeinsame Verteilung der Kurse 
# und des Minimums (S_1,...,S_n, m_n) zu den Zeitpunkten T=(t_1,...,t_n)
# von der normalen GBM, d.h. mit der multivariaten Normalverteilung
# Das Minimum m_n muss dann auch für den Pfadverlauf (S_1,...,S_n) bestimmt werden als m_n|S_1,...,S_n.
# Eine einfache Simulation m_n | S_n reicht nicht!
#
# @param N     Anzahl der Samples/Monte Carlo draws
# @param S
# @param T     Vektor der Bewertungstage T=(t_1,...,t_n) der Länge n
# @param mu    Driftterm
# @param sigma Volatilität
simPricesAndMinimumFromGBM <- function(N=100, S, T, mu, sigma, log=FALSE, m=Inf) {
	
  # Anzahl der Bewertungstage
  n <- length(T)
  
  # Mittelwertvektor mean und Kovarianzmatrix Sigma für die GBM/rmvnorm() bestimmen
  mean  <- (mu - sigma^2/2) * T
  Sigma <- outer(T, T, pmin) * sigma^2
	
  # 1. Simuliere die realisierten Renditen an den Bewertungstagen (r_1,...,r_tk) aus der multivariaten Normalverteilung
  r <- rmvnorm(n=N, mean=mean, sigma=Sigma)
  
  # 2. Simuliere die Periodenminima zwischen den Bewertungstagen aus der Brownschen Brücke
  # Matrix der simulierten Perioden-Minima
  period_min <- matrix(NA, N, n)
  
  for(i in 1:n) {
	  if (i == 1) {
		  period_min[,i] <- rBrownianBridgeMinimum(N, t0=0, T=T[i], a=0, b=r[,i], sigma=sigma)
	  }
	  else {
		  period_min[,i] <- rBrownianBridgeMinimum(N, t0=T[i-1], T=T[i], a=r[,i-1], b=r[,i], sigma=sigma)
	  }
  }
  
  running_min <- matrix(NA, N, n)
  running_min[,1] <- pmin(period_min[,1], m)
  for (i in 2:n) {
    running_min[,i] <- pmin(running_min[,i-1], period_min[,i], m)
  }
  
  # Return matrix ist N x (2*n)  
	R <- cbind(r, running_min)
	if (log == FALSE) {
	  # (S(t_1),...,S(t_n),m(t_1),...,m(t_n))
	  res           <- S * exp(R)
	  colnames(res) <- c(paste("S",1:n, sep=""),paste("m",0:(n-1),1:n,sep=""))
  } else {
	  res           <- R
	  colnames(res) <- c(paste("r",1:n, sep=""),paste("m",0:(n-1),1:n,sep=""))	
	}
  return(res)
}
 

# Simuliere Samples aus der (k+1)-dimensionalen gemeinsame Verteilung der Kurse und des Minimums (S_1,...,S_n, m_n) zu den Zeitpunkten T=(t_1,...,t_n)
# aus einer gestutzten GBM (z.B. S_1 < X_1, S_2 < X_2 usw.) 
# Die Simulation von (S_1,...,S_n) erfolgt mit der gestutzten Normalverteilung.
#
# Das Minimum m_n muss dann auch für den Pfadverlauf (S_1,...,S_n) bestimmt werden als m_n|S_1,...,S_n.
# Eine einfache Simulation m_n | S_n reicht nicht!
#
# Mit Hilfe dieser Funktion kann das Barrierenrisiko von Expresszertifikaten beurteilt werden:
# Bestimme das Barrierenrisiko unter Bedingung, dass das Laufzeitende erreicht wird, d.h. S_1 < X_1, S_2 < X_2, ..., S_{n-1} < X_{n-1}
#
# @param N     Anzahl der Samples
# @param S
# @param T     Vektor der Bewertungstage der Länge k : T=(t_1,...,t_n)'
# @param mu    Driftterm
# @param sigma Volatilität
# @param lowerX Vektor der unteren Kursschwellen der Länge k 
# @param upperX Vektor der oberen Kursschwellen der Länge k
# @param log   
#
# @return entweder 
# log = FALSE : a) die Kurse    (S_1,...,S_n) und die Minimum (m_1,...,m_n)für Expresszertifikate
# log = TRUE  : b) die Renditen (r_1,...,r_n) und das Minimum (m_1,...,m_n) für Expresszertifikate
#
simPricesAndMinimumFromTruncatedGBM <- function(N=100, S, T, mu, sigma, 
  lowerX=rep(0,length(T)), upperX=rep(+Inf, length(T)), log=FALSE, m=Inf) {
	
	# Anzahl der Bewertungstage
	n <- length(T)
	
	# Mittelwertvektor mean und Kovarianzmatrix Sigma für die GBM/rmvnorm() bestimmen
	mean  <- (mu - sigma^2/2) * T
	Sigma <- outer(T, T, pmin) * sigma^2
	
	# 1. Simuliere die realisierten Renditen an den Bewertungstagen (r_1,...,r_tk) aus der gestutzten multivariaten Normalverteilung
	r <- rtmvnorm(n=N, mean=mean, sigma=Sigma, lower=log(lowerX/S), upper=log(upperX/S))
	
	# 2. Simuliere die Periodenminima zwischen den Bewertungstagen aus der Brownschen Brücke
	# Matrix der simulierten Perioden-Minima
	period_min <- matrix(NA, N, n)
	
	for(i in 1:n) {
		if (i == 1) {
			period_min[,i] <- rBrownianBridgeMinimum(N, t0=0, T=T[i], a=0, b=r[,i], sigma=sigma)
		}
		else {
			period_min[,i] <- rBrownianBridgeMinimum(N, t0=T[i-1], T=T[i], a=r[,i-1], b=r[,i], sigma=sigma)
		}
	}
	
	# 3. Bestimme das Gesamtminimum
	#running_min <- apply(m_sim, 1, min)
	running_min <- matrix(NA, N, n)
  running_min[,1] <- pmin(period_min[,1], m)
  for (i in 2:n) {
    running_min[,i] <- pmin(running_min[,i-1], period_min[,i], m)
  }
		
	# Rückgabematrix ist N x (2*n)
	R <- cbind(r, running_min)
	if (log == FALSE) {
		# (S(t_1),...,S(t_n),m(t_1),...,m(t_n))
	  res           <- S * exp(R)
	  colnames(res) <- c(paste("S",1:n, sep=""),paste("m",0:(n-1),1:n,sep=""))  
	}
	else {
		res           <- R
	  colnames(res) <- c(paste("r",1:n, sep=""),paste("m",0:(n-1),1:n,sep=""))
	}
	return(res)
}



# Test der Methode im Vergleich zum Monte-Carlo-Euler-Scheme
if (FALSE) {

  # Simulation der Kurse und des Minimums m_t über direktem Weg	
  X1 <- simulatePricesAndMinimumFromGBM(N=5000, S=100, T=c(1,2,3), mu=0.05, sigma=0.3)
  m1 <- X1[,4]
  
  # Simulation der Kurse und des Minimums über Monte Carlo-Euler-Scheme
  mc.loops <- 5000
  mc.steps <- 2000
  Rprof("test1.out")
  S <- matrix(NA, mc.loops, mc.steps + 1)
  for (j in 1:mc.loops)
  {
	# Simuliere den kompletten Prozess S(t)
	S[j,] <- GBM(S0=100, mu=0.05, sigma=0.3, T=3, N=mc.steps)
  }	
  m2 <- apply(S, 1, min) # zeilenweises (=pfadweises Minimum)
  
  # Vergleich beider Ergebnisse (Dichtefunktion und Verteilungsfunktion des Minimums m_t) und dem theoretischen Ergebnis
  par(mfrow=c(2,2))
  # a) Dichtefunktion des Minimums einer GBM am Laufzeitende vergleichen
  plot(density(m1, to=100), col="black")
  lines(density(m2, to=100), col="blue")
  
  # b) Verteilungsfunktion von m_t vergleichen zwischen beiden Verfahren und dem theoretischen Ergebnis
  B  <- seq(0, 100, by=1)
  p3 <- calculateGBMProbability(Type="P(m_t <= B)", S0=100, B=B, t=3, mu=0.05, sigma=0.3)
  
  plot(ecdf(m1), col="black", main="Vergleich der Minimumverteilung m_t")
  lines(ecdf(m2), col="blue")
  lines(B, p3, col="red")
	legend("topleft", legend=c("Finite-dimensions and Brownian Bridge", "Monte Carlo Euler scheme", "theoretical value"), 
    col=c("black","blue","red"), lwd=2) 
	# --> Passt, in Beispiel übernehmen!
}