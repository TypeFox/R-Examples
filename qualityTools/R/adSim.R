
# Funktion adSim() intialisieren mit Standardeinstellungen: Verteilungsannahme: normal; Bootstrap-Simulationen: 10000

		adSim <- function(x, distribution = "normal", b = 10000){
	
	
# Fehlerpruefungen durchfuehren

# 1.20 ueberpruefung, dass x ein numerischer Vektor ist
		if(mode(x) != "numeric")
			stop(paste("\n","          adSim() requires numeric x data"))	
			
# 1.25 ueberpruefung von x auf fehlende Werte
		if(any(is.na(x))) 
			stop(paste("\n","          x data has missing values (NA)"))
	
# 1.30 ueberpruefung, ob auf die angegebene Verteilung getestet werden kann
		
			# Vektor "distr" enthaelt die testbaren Verteilungsmodelle
			distr <- c("exponential","cauchy","gumbel","gamma","log-normal","lognormal","logistic","normal","weibull")
			
			if(any(distr == distribution)  == FALSE)
				stop(paste("\n","          adSim() can not apply for",distribution,"distribution.				Please choose one of the following distributions for testing goodness-of-fit:				exponential, cauchy, gumbel, gamma, log-normal, lognormal, logistic, normal, weibull"))


# 1.35 & 1.40 ueberpruefung, ob die gewaehlte Verteilung fuer Werte von x<=0 definiert ist bzw. der MLE berechnet werden kann, wenn x entsprechende Werte hat
			
			if(any(x<=0)){
				if(distribution == "exponential" || distribution == "lognormal" || distribution == "log-normal" || distribution == "gamma" || distribution == "weibull"){ 
				 	stop(paste("\n","          adSim() can not apply for", distribution ,"distribution while x contains negative values or x has values equal zero."))}}
				 	
				

# 1.45 & 1.50: bei log-Normalverteilung wird die Stichprobe logarithmiert und anschliessend auf Normalverteilung getestet
# 1.52 "testDistr" gibt die Verteilung an, auf die im Funktionsablauf getestet werden soll (nicht immer identisch mit der hypothetischen Verteilung des Anwenders)
	
		if(distribution == "lognormal" || distribution == "log-normal"){
				x <- log(x)
				testDistr <- "normal"    				}

# 1.54 & 1.56 & 1.58: bei Weibull-Verteilung wird die Stichprobe tranformiert x=-log(x) und anschliessend auf Gumbel-Verteilung getestet
# obwohl die Weibull-Verteilung fuer x>=0 definiert ist, muss aufgrund der Logarithmierung x>0 gelten (die Pruefung erfolgt bereits mit 1.35 & 1.40)

	if(distribution == "weibull"){
				x <- -log(x)
				testDistr <- "gumbel"}

# 1.60 bei den anderen Verteilungen erfolgt der Test an ihrer eigenen Verteilung	
	
	if(distribution != "lognormal" & distribution != "log-normal" & distribution != "weibull"){
				testDistr <- distribution}

# 1.70 Stichprobenumfang bestimmen
	
		n = length(x)
		
# 1.80 Pruefung auf eine Mindeststichprobengroesse
	
		if(n<3)
			stop(paste("\n","          adSim() can not apply for sample sizes n < 3."))

# 1.90 Ordnungsstatistik bilden (aufsteigende Sortierung)
	
		x = sort(x, decreasing = FALSE)	


###################################################################
# 2.00 Parameterschaetzung
###################################################################


# 2.10 & 2.15 Zum Schaetzen der Parameter standardmaessig den MLE benutzen (ausser fuer die Normal-, Cauchy- und Gumbel-Verteilung)
# Parameter der Verteilung aus Stichprobe schaetzen und als Liste "parafit" speichern
	
	if(testDistr != "normal" & testDistr != "gumbel" & testDistr != "cauchy"){
			#library(MASS)
	  
			parafit <- MASS::fitdistr(x,testDistr)
			}

# 2.30 & 2.35 Schaetzen der Parameter der Normalverteilung 
# der MLE ist fuer die Standardabweichung der Normalverteilung nicht erwartungstreu
# Parameter der Normalverteilung werden mit mean() und sd()  geschaetzt und als Vektor "parafit" gespeichert	

	if(testDistr == "normal"){
				parafit 	<- numeric(2)
				parafit[1] 	= mean(x) 
				parafit[2] 	= sd(x)
				
				# 2.40 Parameterbenennung
				if(distribution == "lognormal" || distribution == "log-normal"){
					names(parafit) = c(  "meanlog",  "sdlog")   # 2.41
    			}else{names(parafit) = c(  "mean",  "sd")}		# 2.42
			}
	
# 2.50 Schaetzen der Cauchy-Parameter
# Fuer die Cauchy-Verteilung sind zwei unterschiedliche Schaetzverfahren implementiert

	if(testDistr == "cauchy"){

		# 2.52 Pruefung, ob simuliert wird oder kritische Werte ausgelesen werden sollen	
		
		if(is.na(b) == FALSE){		# die AD-Verteilung soll spaeter simuliert werden
			
			# 2.55 Schaetzung mittels fitdistr()	
				
				#library(MASS)
				parafit <- MASS::fitdistr(x,testDistr)
				
		}else{						# 2.52 die AD-Verteilung soll spaeter nicht simuliert werden			
			
		# 2.60 Parameterschaetzung, basierend auf den Summen der gewichteten Ordnungsstatistik
						
			parafit 	<- numeric(2)
		
			uWeight = numeric(n)
			uWeight[1:n] <- sin( 4*pi* ( 1:n/(n+1) - 0.5)  ) / (n*tan( pi* ( 1:n/(n+1) - 0.5 )  ) )  # Berechnung der Gewichte u
		
			#  fuer ungerade n koennen nicht alle Gewichte berechnet werden, deshalb wird das fehlende Gewicht geschaetzt
			if(ifelse(n %% 2, TRUE, FALSE)){ 	# ifelse-Bedingung gibt TRUE, wenn n ungerade ist
				if(length(na.omit(uWeight)) + 1 == length(uWeight)){  # Bedingung prueft, ob nur ein Wert in uWeight nicht berechnet werden konnte
				
				 uWeight[which(is.na(uWeight))] = (n+1)/n - sum(na.omit(uWeight))   # Ersetzen von NaN durch geeigneten Schaetzer 
				}}				
			parafit[1] 	<- uWeight %*% x 			# Parameterberechnung

			
			vWeight = numeric(n)
			vWeight[1:n] <- 8*tan( pi* ( 1:n/(n+1) - 0.5)  )*(cos( pi*(1:n/(n+1) - 0.5 )))^4 / n  # Berechnung der Gewichte v
			
			parafit[2]	<- vWeight %*% x 			# Parameterberechnung		
			
			# 2.65 Parameterbenennung
			names(parafit) = c(  "location",  "scale")
			}
			}
	

# 2.70 Schaetzen der Gumbel-Parameter (fitdistr() kann nicht die Gumbel-Verteilung schaetzen)

	if(testDistr == "gumbel"){
		
				### 2.72 p-Gumbel-Funktion definieren (aus Library VGAM uebernommen)
						
					pgumbel = function (q, location = 0, scale = 1){
    								answer = exp(-exp(-(q - location)/scale))
   									answer[scale <= 0] = NaN
    								answer}
				
				### 2.75 r-Gumbel-Funktion definieren (aus Library VGAM uebernommen)
					
					is.Numeric = function (x, allowable.length = Inf, integer.valued = FALSE, positive = FALSE){
									if (all(is.numeric(x)) && all(is.finite(x)) && (if (is.finite(allowable.length)) length(x) == 
    									allowable.length else TRUE) && (if (integer.valued) all(x == round(x)) else TRUE) &&
    									(if (positive) all(x > 0) else TRUE)) TRUE else FALSE}

					rgumbel = function (n, location = 0, scale = 1){
    								use.n = if ((length.n <- length(n)) > 1)
        							length.n
   									else if (!is.Numeric(n, integ = TRUE, allow = 1, posit = TRUE))
        							stop("bad input for argument 'n'")
    								else n
    								answer = location - scale * log(-log(runif(use.n)))
    								answer[scale <= 0] = NaN
    								answer}
			
				
				### 2.80 MLE-Parameterschaetzung fuer die Gumbel-Verteilung (iterative Loesung)
			
				f <- function(p) 1/n*sum(x) - (sum(x*exp(-x/p)))/(sum(exp(-x/p)))-p
					
					### iterative Loesung notwendig, um beta zu schaetzen (D'Agostino/Stephens 1986, S. 146)
				
					itSol <- uniroot(f,c(-100,100),tol = 0.0000001, maxiter = 100000)
					beta <- as.numeric(itSol$root)
					alpha = - beta * log(sum(exp(-x/beta)/n))

				parafit 	<- numeric(2)
				parafit[1] 	= alpha 
				parafit[2] 	= beta
				
				# 2.85 Parameterbenennung
				names(parafit) = c(  "location",  "scale")
							
			}
	


###################################################################
#          3.00 AD-Teststatistik berechnen (mit Lewis-Formel)
###################################################################

	# 3.10 Transformation des Wahrscheinlichkeitsintegrals (PIT) berechnen

		pit = numeric(n)
	
			#3.15		
			if(testDistr == "normal")
					pit[1:n] = pnorm(x[1:n],parafit[1],parafit[2])
 					
			#3.20
			if(testDistr == "exponential")
 					pit[1:n] = pexp(x[1:n],parafit$estimate["rate"])

			#3.25 							
	 		if(testDistr == "cauchy")
				if(is.na(b) == FALSE){
						pit[1:n] = pcauchy(x[1:n],parafit$estimate["location"],parafit$estimate["scale"])
				}else{  pit[1:n] = pcauchy(x[1:n],parafit[1],parafit[2])}
 					
			#3.30			
			if(testDistr == "gamma")
					pit[1:n] = pgamma(x[1:n],parafit$estimate["shape"],parafit$estimate["rate"])
 					
			#3.35 			
	 		if(testDistr == "gumbel")
					pit[1:n] = pgumbel(x[1:n],parafit[1],parafit[2])

			#3.40	 					
	 		if(testDistr == "logistic")
					pit[1:n] = plogis(x[1:n],parafit$estimate["location"],parafit$estimate["scale"])
 						
 			
 	# 3.50 Anwendung der AD-Formel

			# 3.60 fuer jedes i den Summenwert berechnen und in Matrix h schreiben
			h = matrix(ncol=n)

			h[1,] = (2*col(h)-1)*log(pit) + (2*n + 1 - 2*col(h))* log(1 - pit) 
			
			# 3.70 AD-Formel vollstaendig berechnen
			AD = -n-(1/n) * rowSums(h)
		

			# 3.80 Fehlerpruefung, ob AD-Wert vorhanden ist
			if(is.na(AD) == TRUE || AD<0)
				stop(paste("\n","        The calculation of the Anderson Darling statistic fails."))
	
			# Ende der AD-Wert-Berechnung fuer die beobachteten Daten
			
if(is.na(b) == FALSE){

	
#####################################################################################
#       4 parametrisches Bootstrapping der AD-Nullverteilung
#####################################################################################
	
	# 4.00 Fehlerpruefung, ob die Simulationshaeufigkeit zulaessig gewaehlt ist
	
		if(b<1000)
			stop("b is chosen too small for generate an accurate p-Value.")
		if(b>1000000)
			stop("b is chosen too big for generate an p-value within a reasonable time.")
			
	# 4.02 Ausgabe, dass simuliert wird
		cat("\n","   ... simulating the Anderson-Darling distribution by",b,"bootstraps for",distribution,"distribution...","\n","\n")


	# 4.05 Resampling und Parameterschaetzung
	# Matrix Y mit b sortierten Zufallsstichproben mit Umfang n erstellen
	# Verteilungsparameter zu jeder Bootstrap-Stichprobe schaetzen	
	
								
			# 4.10
			if(testDistr == "normal"){
				Y = t(replicate(b, sort(rnorm(n,parafit[1],parafit[2]))))
				
				paraMean = rowMeans(Y)
				paraSd = sqrt(rowSums((Y-rowMeans(Y))^2) /(ncol(Y)-1))
				}

			# 4.15								
			if(testDistr == "exponential"){
				Y = t(replicate(b, sort(rexp(n, parafit$estimate["rate"]))))
								
 				unParaList <- unlist(apply(Y,1,function(x) fitdistr(x,testDistr,rate = parafit$estimate["rate"])))
				paraRate = unParaList[names(unParaList) == "estimate.rate"]
					}
 			
 			# 4.20			
			if(testDistr == "cauchy"){
				Y = t(replicate(b, sort(rcauchy(n, parafit$estimate["location"],parafit$estimate["scale"]))))
								
				unParaList <- unlist(apply(Y,1,function(x) fitdistr(x,testDistr,list(location = parafit$estimate["location"], scale = parafit$estimate["scale"]))))
				paraLocation = unParaList[names(unParaList) == "estimate.location"]
				paraScale = unParaList[names(unParaList) == "estimate.scale"]
					}
			
			# 4.25
			if(testDistr == "gamma"){
				Y = t(replicate(b, sort(rgamma(n, parafit$estimate["shape"],parafit$estimate["rate"]))))
				
				unParaList <- unlist(apply(Y,1,function(x) fitdistr(x,testDistr,list(shape = parafit$estimate["shape"], rate = parafit$estimate["rate"]))))
				paraShape = unParaList[names(unParaList) == "estimate.shape"]
				paraRate = unParaList[names(unParaList) == "estimate.rate"]
					}
			# 4.30	
			if(testDistr == "gumbel"){
				Y = t(replicate(b, sort(rgumbel(n, parafit[1], parafit[2]))))
								 				
 				paraBeta	= numeric(b)
 				paraAlpha	= numeric(b)
 					
 				for(j in 1:b){
 					 					
					### iterative Loesung notwendig, um beta zu schaetzen (D'Agostino/Stephens 1986, S. 146)
				
					itSol <- uniroot(function(p) 1/n*sum(Y[j,])-
									(sum(Y[j,]*exp(-Y[j,]/p)))/(sum(exp(-Y[j,]/p)))-p,c(-100,100),tol = 0.0000000001, maxiter = 100000)
				paraBeta[j] <- as.numeric(itSol$root)
				paraAlpha[j] = - paraBeta[j] * log(sum(exp(-Y[j,]/paraBeta[j])/n))
						}
				}				
			
			#4.35 			
	 		if(testDistr == "logistic"){
				Y = t(replicate(b, sort(rlogis(n, parafit$estimate["location"],parafit$estimate["scale"]))))
										
				unParaList <- unlist(apply(Y,1,function(x) fitdistr(x,testDistr,list(location = parafit$estimate["location"], scale = parafit$estimate["scale"]))))
				paraLocation = unParaList[names(unParaList) == "estimate.location"]
				paraScale = unParaList[names(unParaList) == "estimate.scale"]
				}
							 			

######################################################################
### 4.40 b-fache Berechnung der Anderson-Darling-Formel
######################################################################

# PIT fuer jeden Wert der Matrix Y berechnen
# (abhaengig von den geschaetzten Verteilungsparameter jeder Bootstrap-Stichprobe)

			#4.45
			if(testDistr == "normal")
					Y[,1:n] <- pnorm(Y[,1:n],paraMean,paraSd)
				
			#4.50			
			if(testDistr == "exponential")
					Y[,1:n] <- pexp(Y[,1:n],paraRate) 		
 		
 			#4.55
 			if(testDistr == "cauchy")
					Y[,1:n] <- pcauchy(Y[,1:n],paraLocation,paraScale)			
 			#4.60
			if(testDistr == "gamma")
					Y[,1:n] <- pgamma(Y[,1:n],paraShape,paraRate)
 		
 			#4.65
 			if(testDistr == "gumbel")
				Y[,1:n] <- pgumbel(Y[,1:n],paraAlpha,paraBeta)

			#4.70
	 		if(testDistr == "logistic")
					Y[,1:n] <- plogis(Y[,1:n],paraLocation,paraScale)


# 4.75 Berechnung der Summenglieder der AD-(Lewis)-Formel zu jeder Bootstrap-Stichprobe (Matrix Y ueberschreiben)

Y[1:b,] <- (2*col(Y)-1)*log(Y[1:b,]) +  (2*n + 1 - 2*col(Y)) * log(1 - Y[1:b,])
 	
# 4.77 simulierte AD-Werte berechnen und in Vektor simAD schreiben

	d = rowSums(Y)
	
	simAD = numeric(b)    

	simAD[1:b] = -n-(1/n)*d[1:b]

# 4.80 Fehlerpruefung, ob alle AD-Werte bestimmt werden konnten
if(any(is.na(simAD))){
	cat("    The simulated Anderson-Darling distribution contains NAs or NaNs!","\n","\n")}
		

# 4.90 Bestimmung kritischer Werte

critValues = round (matrix(	c(  0.75, 0.90, 0.95, 0.975, 0.990, 
				quantile(simAD, 0.75,na.rm = TRUE), quantile(simAD, 0.90,na.rm = TRUE), quantile(simAD, 0.95,na.rm = TRUE), 
				quantile(simAD, 0.975,na.rm = TRUE), quantile(simAD, 0.99,na.rm = TRUE) ), nrow= 2, byrow = TRUE ), digits = 5)


# 4.95 Bestimmung des p-Werts
# die p-Wert Bestimmung ist so programmiert, dass fehlende Werte in simAD ausgelassen werden (na.omit)

pValue = sum(na.omit(simAD) > AD)/length(na.omit(simAD))

# Ende der simulationsbasierten Bestimmung von kritischen Werten bzw. des p-Werts


}else{
	
		# 5.00 (b = NA), d.h. Auslesen tabellierter kritischer Werte bzw. des p-Werts oder Anwendung von Berechnungsformeln zur p-Wert-Bestimmung
		
		# 5.02 simAD wird spaeter in der Ausgabe abgerufen, weshalb die Variable definiert sein muss
		simAD = NA
		
		### 5.06 Definieren einer Matrix "critValues" zur Erfassung kritischer Werte
		
		critValues = matrix(	c(  0.75, 0.90, 0.95, 0.975, 0.990, NA, NA ,NA, NA, NA), nrow = 2, byrow = TRUE)
			
		
		# 5.08 Test auf Normal-, Exponential-, Gumbel- oder logistische Verteilung
			if(testDistr == "normal" || testDistr == "exponential" ||testDistr == "gumbel" || testDistr == "logistic"){
					
			# 5.10 Test auf Normalverteilung
				if(testDistr == "normal"){
						
					###################################################################################################################
					### kritische Werte fuer die Normalverteilung nach D'Agostino und Stephens 1986 - S. 123   Tab. 4.7 
					################################################################################################################### 					
					
					# 5.12 Tabelle der kritischen Werte
					normalMtx = matrix(c(	0.50 , 0.75, 0.85 , 0.90 , 0.95 , 0.975 , 0.99 , 0.995 , 
											0.341 , 0.470 , 0.561 , 0.631 , 0.752 , 0.873 , 1.035 , 1.159), nrow=2, byrow = TRUE )
					
					# 5.14 Anpassen der Wertetabelle bezueglich des Stichprobenumfangs
					normalMtx[2,1:ncol(normalMtx)] = normalMtx[2,1:ncol(normalMtx)]/(1+0.75/n+2.25/n^2)
					refMtx = normalMtx
					
					########################################################################################################
					### p-Wert Berechnung fuer die Normalverteilung nach D'Agostino und Stephens 1986 - S. 127    
					######################################################################################################## 
					
					# 5.16 Anpassen des AD-Werts bzgl. des Stichprobenumfangs
					cAD = AD*(1+0.75/n+2.25/n^2)
					
					# 5.18 Berechnung des p-Werts nach Formel
					if(0.600<cAD)
						{pValue = exp(1.2937-5.709*cAD+0.0186*cAD^2)}
					if(0.340<cAD & cAD<0.600)
						{pValue = exp(0.9177-4.279*cAD-1.38*cAD^2)}
					if(0.200<cAD & cAD<0.340)
	 					{pValue = 1-exp(-8.318+42.796*cAD-59.938*cAD^2)}
					if(cAD<0.200)
						{pValue = 1-exp(-13.436+101.14*cAD-223.73*cAD^2)}
				}	
					
				# 5.20 Test auf Exponentialverteilung
				if(testDistr == "exponential"){
					
					#############################################################################################################
					### Kritische Werte fuer die Exponentialverteilung nach D'Agostino und Stephens 1986 - S. 135   Tab. 4.11   
					############################################################################################################# 					
					
					# 5.22 Tabelle der kritischen Werte
						expMtx = matrix( c( 0.75, 0.80, 0.85 , 0.90 , 0.95 , 0.975 , 0.99 , 0.995 , 0.9975,
										0.736 , 0.816 , 0.916 , 1.062 , 1.321 , 1.591 , 1.959 , 2.244, 2.534), nrow = 2, byrow = TRUE )
					
					# 5.24 Anpassen der Wertetabelle bezueglich des Stichprobenumfangs
						expMtx[2,1:ncol(expMtx)] = expMtx[2,1:ncol(expMtx)]/(1+0.6/n)
	
						refMtx = expMtx
				
					############################################################################################################
					### p-Wert-Berechnung fuer die Exponentialverteilung nach D'Agostino und Stephens 1986 - S. 136
					### Case 2: origin (bzw. location) = 0 = known , scale = unknown = 1/Lambda = 1/rate  
					############################################################################################################
					
					# 5.26 Anpassen des AD-Werts bzgl. des Stichprobenumfangs
						cAD = AD*(1+0.6/n)
					
					# 5.28 Berechnung des p-Werts nach Formel
					if(0.950<cAD)
						{pValue = exp(0.731 - 3.009*cAD + 0.15*cAD^2)}
					if(0.510<cAD & cAD<0.950)
						{pValue = exp(0.9209 - 3.353*cAD + 0.300*cAD^2)}
					if(0.260<cAD & cAD<0.510)
						{pValue = 1 - exp(-6.1327 + 20.218*cAD - 18.663*cAD^2)}
					if(cAD<0.260)
						{pValue = 1 - exp(-12.2204 + 67.459*cAD - 110.3*cAD^2)}				
				}		
				
			# 5.30 Test auf logistische Verteilung oder Gumbel-Verteilung
			if(testDistr == "gumbel" || testDistr == "logistic"){			
			
				# 5.31 Test auf Gumbel-Verteilung
				if(testDistr == "gumbel"){
					
					##############################################################################################################
					### Kritische Werte fuer die Extremwert-Verteilung nach D'Agostino und Stephens 1986 - S. 146   Tab. 4.17   
					############################################################################################################## 					
					
					# 5.32 Tabelle der kritischen Werte
						gumbelMtx = matrix(c( 0.75, 0.90 , 0.95 , 0.975 , 0.99 , 0.474, 0.637, 0.757, 0.877, 1.038), nrow = 2, byrow = TRUE )
					
											
					# 5.34 Anpassen der Wertetabelle bezueglich des Stichprobenumfangs
						gumbelMtx[2,1:ncol(gumbelMtx)] = gumbelMtx[2,1:ncol(gumbelMtx)]/(1 + 0.2/sqrt(n))
						refMtx = gumbelMtx
				}
			
				# 5.35 Test auf logistische Verteilung			
				if(testDistr == "logistic"){
						
					##############################################################################################################
					### kritische Werte fuer die logistische Verteilung nach D'Agostino und Stephens 1986 - S. 157   Tab 4.22 
					############################################################################################################## 
										
					# 5.37 Tabelle der kritischen Werte
						logisMtx = matrix(c( 0.75, 0.90 , 0.95 , 0.975 , 0.99 , 0.995 , 0.426, 0.563, 0.660, 0.769, 0.906, 1.010 ), nrow = 2, byrow = TRUE )
					
					# 5.39 Anpassen der Wertetabelle bezueglich des Stichprobenumfangs
						logisMtx[2,1:ncol(logisMtx)] = logisMtx[2,1:ncol(logisMtx)]/(1 + 0.25/n)
						refMtx = logisMtx
		}
			
						
				##############################################################################################################
				### Bestimmung des p-Werts fuer die Gumbel- oder logistische Verteilung
				############################################################################################################## 
			
					critCheck <- refMtx[2,1:ncol(refMtx)] > AD   ### gibt TRUE aus fuer alle Eintraege in der Zeile, die groesser als AD sind
					
						# 5.40 Existiert ein kritischer Wert, der groesser als der AD-Wert ist?   	
						if(any(critCheck)){
							
							
							# 5.42 firPos gibt die Position in der Zeile vom letzten Wert an, der noch kleiner als der AD-Wert ist
								firPos <- min(which(critCheck)) - 1

						}else{
							# 5.44 letzte Spalte als Position waehlen
								firPos <- ncol(refMtx)}
							
							# 5.46 p-Wert entsprechend der ermittelten Position bestimmen
								if(firPos == 0){
									pValue <- 1 - refMtx[1,1]
									pValue <- paste(">",pValue) 
								}else{
									pValue <- 1 - refMtx[1,firPos]
									pValue <- paste("<=",pValue)}}
										
				##############################################################################################################
				### Auslesen der kritischen Werte fuer Normal-, Exponential-, Gumbel oder logistische Verteilung				
				############################################################################################################## 
												
					for(i in 1:ncol(critValues)){
						
						# 5.50 Ist der kritische Werte fuer das zu spezifizierende Quantil tabelliert?
						if(any(refMtx[1,1:ncol(refMtx)] == critValues[1,i])){                          
							
								
							# 5.52 dann Position des kritischen Werts bzgl. des betrachteten Quantils erfassen
								position <- (1:length(refMtx[1,1:ncol(refMtx)] == critValues[1,i]))[(refMtx[1,1:ncol(refMtx)] == critValues[1,i])] 
							
							# 5.54 Auslesen des kritischen Werts bezueglich der Position	
								critValues[2,i] <- refMtx[2,position]          ### liest aus der Matrix den kritschen Wert abhaengig vom gewaehlten Quantil
						}else{
							# 5.58 nicht-tabellierte Quantile mit "NA" belegen
								critValues[2,i] <- NA}}
				}			
			
	# 5.60 Test auf Gamma-Verteilung
	if(testDistr == "gamma"){
					
					###################################################################################################################################
					### Bestimmung des kritischen Werts und des p-Werts fuer die Gammaverteilung nach D'Agostino und Stephens 1986 - S. 155 - Tab. 4.21
					### case 3: shape = unknown, scale = unknown = 1/rate, origin = known = 0    
					###################################################################################################################################
					
					# 5.62 Tabelle der kritischen Werte
					gammaDF = data.frame(	c( 1,     2,     3,     4,     5,     6,     8,    10,    12,    15,    20,     Inf  ),
											c( 0.486, 0.477, 0.475, 0.473, 0.472, 0.472, 0.471, 0.471, 0.471, 0.47,  0.47,  0.47 ), 
											c( 0.657, 0.643, 0.639, 0.637, 0.635, 0.635, 0.634, 0.633, 0.633, 0.632, 0.632, 0.631),
						 			 	 	c( 0.786, 0.768, 0.762, 0.759, 0.758, 0.757, 0.755, 0.754, 0.754, 0.754, 0.753, 0.752),
											c( 0.917, 0.894, 0.886, 0.883, 0.881, 0.88,  0.878, 0.877, 0.876, 0.876, 0.875, 0.873),
											c( 1.092, 1.062, 1.052, 1.048, 1.045, 1.043, 1.041, 1.04,  1.039, 1.038, 1.037, 1.035),
						 					c( 1.227, 1.19,  1.178, 1.173, 1.17,  1.168, 1.165, 1.164, 1.163, 1.162, 1.161, 1.159))

					names(gammaDF) = c(  "m",  0.75,  0.90,  0.95,  0.975, 0.99,  0.995)

							
					###################################################################################################################################
					######## p-Wert-Bestimmung fuer die Gamma-Verteilung
					###################################################################################################################################
					
						# critCheck gibt TRUE aus fuer alle Eintraege in der entsprechenden m-Zeile, die groesser als der AD-Wert sind
						# zu beachten ist, dass die Betrachtung erst ab der 2ten Spalte der Tabelle erfolgt, so dass die Indizierung von critCheck
						# gegenueber der Indizierung der Tabelle  um -1 verschoben ist 
						
							critCheck <- gammaDF[min(which(gammaDF$m >= parafit$estimate["shape"])),2:ncol(gammaDF)] > AD   
								 
						# 5.65 Existiert ein kritischer Wert, der groesser als der AD-Wert ist (in der entsprechenden m-Zeile)?  
							if(any(critCheck)){	
								
								# 5.66 firPos gibt die Spalten-Position der Tabelle in der entsprechenden m-Zeile von dem krit. Wert an, 
								# der noch kleiner als der AD-Wert ist. Fuer firPos = 1 existiert kein kleiner kritischer Wert
								
									firPos <- min(which(critCheck))     
							}else{
								
							# 5.67 letzte Spalte als Position waehlen
									firPos <- ncol(gammaDF) }
							
							# 5.68  p-Wert entsprechend der ermittelten Position bestimmen	
							if(firPos == 1){							
									pValue <- 1 - as.numeric(names(gammaDF)[2])       
									pValue <- paste(">",pValue) 
								}else{
									pValue <- 1 - as.numeric(names(gammaDF)[firPos])     
									pValue <- paste("<=",pValue)}

					###################################################################################################################################
					######## kritischen Wert fuer die Gamma-Verteilung auslesen
					###################################################################################################################################
						  
							
							for(i in 1:ncol(critValues)){
								
								# 5.70 Ist der kritische Wert fuer das zu spezifizierende Quantil tabelliert?
								if(any(names(gammaDF) == critValues[1,i] )){
									
									
									# 5.72 Auslesen der kritischen Werte an der Zeilen-Position des entsprechenden Formparameters
									# und der Spalten-Position des zu spezifizierenden Quantils
										critValues[2,i] <- gammaDF[min(which(gammaDF$m >= parafit$estimate["shape"] )),which(names(gammaDF) == critValues[1,i])]
								}else{
									# 5.74 nicht-tabellierte Quantile mit "NA" belegen
										critValues[2,i] <- NA}}
					}
	
	# 5.80 Test auf Cauchy-Verteilung			
	if(testDistr == "cauchy"){
					
					####################################################################################################################
					### Bestimmung des kritischen Werts fuer die Cauchy-Verteilung nach D'Agostino und Stephens 1986 - S. 163   Tab 4.26
					### case 3: location = unknown, shape = unknown        
					#################################################################################################################### 					
					
					# 5.82 Tabelle der kritischen Werte
							cauchyDF = data.frame(
									c( 5,     8,     10,   12,   15,   20,    25,    30,    40,    50,    60,    100,   Inf),
									c( 0.835, 0.992, 1.04, 1.04, 1.02, 0.975, 0.914, 0.875, 0.812, 0.774, 0.743, 0.689, 0.615),
									c( 1.14,  1.52,  1.63, 1.65, 1.61, 1.51,  1.4,   1.3,   1.16,  1.08,  1.02,  0.927, 0.78),
									c( 1.4,   2.06,  2.27, 2.33, 2.28, 2.13,  1.94,  1.76,  1.53,  1.41,  1.3,   1.14,  0.949),
									c( 1.77,  3.2,   3.77, 4.14, 4.25, 4.05,  3.57,  3.09,  2.48,  2.14,  1.92,  1.52,  1.225),
									c( 2,     4.27,  5.58, 6.43, 7.2,  7.58,  6.91,  5.86,  4.23,  3.37,  2.76,  2.05,  1.52),
									c( 2.16,  5.24,  7.5,  9.51, 11.5, 14.57, 14.96, 13.8,  10.2,  7.49,  5.32,  3.3,   1.9))
											
					names(cauchyDF) = c(  "n",  0.75,  0.85, 0.90,  0.95,  0.975, 0.99)						

								
					##############################################################################################################
					### Bestimmung des p-Werts fuer die Cauchy-Verteilung 
					############################################################################################################## 
				
					# 5.84 Existieren kritische Werte fuer den zu untersuchenden Stichprobenumfang n?
					if(any(cauchyDF[1:13,1] == n)){						
					
						# critCheck gibt TRUE aus fuer alle Eintraege in der entsprechenden n-Zeile, die groesser als der AD-Wert sind
						# zu beachten ist, dass die Betrachtung erst ab der 2ten Spalte der Tabelle erfolgt, so dass die Indizierung von critCheck
						# gegenueber der Indizierung der Tabelle um -1 verschoben ist				
								critCheck <- cauchyDF[which(cauchyDF[1:13,1] == n),2:ncol(cauchyDF)] > AD
						
						# 5.85 Existiert ein kritischer Wert, der groesser als der AD-Wert ist (in der entsprechenden n-Zeile)? 
							if(any(critCheck)){	
									
								# 5.86 firPos gibt die Spalten-Position der Tabelle in der entsprechenden n-Zeile von dem krit. Wert an, 
								# der noch kleiner als der AD-Wert ist; fuer firPos = 1 existiert kein kleiner kritischer Wert
									firPos <- min(which(critCheck))
								
							}else{
								# 5.87  letzte Spalte als Position waehlen 	
									firPos <- ncol(cauchyDF) }
							
						# 5.88  p-Wert entsprechend der ermittelten Position bestimmen
							if(firPos == 1){
									pValue <- 1 - as.numeric(names(cauchyDF)[2])
									pValue <- paste(">",pValue) 
								}else{
									pValue <- 1 - as.numeric(names(cauchyDF)[firPos])
									pValue <- paste("<=",pValue)}
									
					##############################################################################################################
					### Kritische Werte fuer die Cauchy-Verteilung auslesen 
					############################################################################################################## 
					
											
							for(i in 1:ncol(critValues)){
								
								# 5.90 Ist der kritische Wert fuer das zu spezifizierende Quantil tabelliert?
								if(any(names(cauchyDF) == critValues[1,i] )){
																		
									# 5.92 Auslesen der kritischen Werte an der Zeilen-Position des entsprechenden Stichprobenumfangs n
									# und der Spalten-Position des zu spezifizierenden Quantils
										critValues[2,i] <- cauchyDF[which(cauchyDF[1:13,1] == n),which(names(cauchyDF) == critValues[1,i])]
								}else{
									# 5.94 nicht-tabellierte Quantile mit "NA" belegen
										critValues[2,i] <- NA}
							}			
																		
												
					}else{						
						# 5.96 p-Wert kann nicht ermittelt werden, da fuer "n" keine tabellierten Werte existieren
							pValue <- NA
						
						# 5.98 kritische Werte koennen nicht ermittlet werden, da fuer "n" keine tabellierten Werte existieren
							critValues[2,1:ncol(critValues)] <- NA
							cat("\n","Critical values / p-Values for the Cauchy Distribution are only tabled for sample sizes: n = 5, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 100","\n")}
											
		}								
}

###########################################################################
####	6.00 Parameterkorrektur fuer die Weibull-Verteilung
###########################################################################

# 6.10 Sind die beobachteten Zufallszahlen auf Weibull-Verteilung getestet worden?
if(distribution == "weibull"){
	
	# 6.20 Umbenennung der Parameternamen entsprechend der Weibull-Verteilung
		names(parafit) = c(  "shape",  "scale")

	# 6.30 Kopie des "parafit"-Vektors 
		parafitCopy = parafit
			
	# 6.40 Umrechnung der Gumbel-Parameter in Weibull-Parameter
    	parafit[1] 	= (1/parafitCopy[2])        
		parafit[2] 	= exp(- parafitCopy[1])}


print(list(distribution = distribution, parameter_estimation = parafit,Anderson_Darling = AD, p_value = pValue))
invisible(list(distribution = distribution, parameter_estimation = parafit,Anderson_Darling = AD,p_value = pValue,crititical_values = critValues,simAD = simAD))

}


		