##############################################################################################
#                                                                           
# PROCEDURA DI REGIONALIZZAZIONE DELLE PORTATE DI PIENA ARPIEM 2012
#                                                                           
# AUTHOR(S):    Ganora Daniele                                           
#                                                                         
# PURPOSE:      Raccolta di funzioni R per l'applicazione della procedura ARPIEM 2012.
#			             
# VERSION:		jun-2013
#
##############################################################################################

##############################################################################################
# Stima degli L-momenti campionari
ARPIEM2012.Lmom.sample <- function(x){
	# caricamento campione (solo dati sistematici)
	x <- as.matrix(x)

	# calcolo Lmomenti e relative deviazioni standard
	n <- length(x)
	lmoms <- Lmoments(x)
	Qind <- lmoms["l1"]
	LCV <- lmoms["lcv"]
	LCA <- lmoms["lca"]
	# Lkur.C <- lmoms["lkur"]

	# deviazione standard Qind
	sigma.Qind <- sqrt( var(x) / n )
	
	# deviazione standard secondo Elamir & Seheult (2004)
	sigma.LCV.ES <- sqrt( nsRFA::varLCV(x) )
	sigma.LCA.ES <- sqrt( nsRFA::varLCA(x) )
	
	# deviazione standard formule semplificate (Viglione, 2007)
	sigma.LCV.V <- sqrt( (0.9 * LCV)^2 / n )
	sigma.LCA.V <- sqrt( (0.45 + 0.6*abs(LCA))^2 / n )

	# dove la formulazione di Elamir & Seheult (2004) fallisce, viene utilizzata quella semplificata di Viglione (2007)	
	if (is.na(sigma.LCV.ES)) {
		sigma.LCV <- sigma.LCV.V
	} else {
		sigma.LCV <- sigma.LCV.ES	
	}

	if (is.na(sigma.LCA.ES)) {
		sigma.LCA <- sigma.LCA.V
	} else {
		sigma.LCA <- sigma.LCA.ES	
	}

	# restituzione valori
	note <- c("C", "C", "C")
	variabile <- c("Qind", "LCV", "LCA")
	media <- c(Qind, LCV, LCA)
	sigma <- c(sigma.Qind, sigma.LCV, sigma.LCA)
	output <- data.frame(variabile, media, sigma, note)
	rownames(output) <- variabile
	return(output)
}
##############################################################################################

##############################################################################################
# Stima degli L-momenti regionali
ARPIEM2012.Lmom.reg <- function(descr){
	options(stringsAsFactors=FALSE)
		
	nomi <- c("area_bacinokm", "quota_minima", "IDFa", "IDFa_cv", "IDFn", "LCV1h", "NDVIanno", "LCV6h_cv", "LCA6h", "LCA24h_cv", "fourier_B2", "clc2" )
	
	controllo <- descr[nomi]; controllo <- controllo[!is.na(controllo)]; if (length(controllo)!=12){print("Errore: ")} else {
		
 	descrittori <- as.matrix(descr, ncol=1)
 	descrittori.log <- log( descrittori[c("area_bacinokm", "IDFa", "IDFn", "LCV1h", "quota_minima", "NDVIanno", "IDFa_cv", "LCV6h_cv"), , drop=F] )
	rownames(descrittori.log) <- paste("log.", rownames(descrittori.log), sep="")
	descrittori <- rbind(intercept=1, descrittori, descrittori.log)
 
	# calcolo valori regionali
	ARPIEM2012regressions <- NULL
	load(system.file("extdata","ARPIEM2012regressions.RData",package="hydroApps"))
		# Qind.R2 ###########################################################################
		reg <- ARPIEM2012regressions$Qind
		xT <- descrittori[reg$model, , drop=F]
				
		Qind.R2 <- exp( t(xT) %*% reg$coefficients )
		var.logQindR2 <- reg$varmodel + t(xT) %*% solve(t(reg$X) %*% solve(reg$Lambda) %*% reg$X) %*% xT
		sigma.Qind.R2 <- sqrt( Qind.R2^2 * (exp(var.logQindR2) - 1) )
		
		# LCV.R1 ###########################################################################
		reg <- ARPIEM2012regressions$LCV
		xT <- descrittori[reg$model, , drop=F]
				
		LCV.R1 <- exp( t(xT) %*% reg$coefficients )
		var.logLCVR1 <- reg$varmodel + t(xT) %*% solve(t(reg$X) %*% solve(reg$Lambda) %*% reg$X) %*% xT
		sigma.LCV.R1 <- sqrt( LCV.R1^2 * (exp(var.logLCVR1) - 1) )
		
		# LCA.R ############################################################################
		reg <- ARPIEM2012regressions$LCA
		xT <- descrittori[reg$model, , drop=F]
				
		LCA.R <-  t(xT) %*% reg$coefficients 
		sigma.LCA.R <- sqrt( reg$varmodel + t(xT) %*% solve(t(reg$X) %*% solve(reg$Lambda) %*% reg$X) %*% xT )

	# restituzione valori
	note <- c("R", "R", "R")
	variabile <- c("Qind", "LCV", "LCA")
	media <- c(Qind.R2, LCV.R1, LCA.R)
	sigma <- c(sigma.Qind.R2, sigma.LCV.R1, sigma.LCA.R)
	output <- data.frame(variabile, media, sigma, note)
	rownames(output) <- variabile
	return(output)
	}
}
##############################################################################################


ARPIEM2012.sim.Lmoments <- function(Qind.type, LCV_LCA.type, Qind, sdQind, LCV, sdLCV, LCA, sdLCA, n=10000){
	# tipo.Qind:
	#	"C" = Qind campionaria
	#	"R" = Qind regionale modello moltiplicativo
	# tipo.LCV_LCA:
	#	"C_C" = LCV campionario e LCA campionario
	#	"C_R" = LCV campionario e LCA regionale modello additivo
	#	"R_R" = LCV regionale modello moltiplicativo e LCA regionale modello additivo 
	
		# Simulazione Qind campionaria
	ARPIEM2012.sim.Qind_camp <- function(mQ, sdQ, n){
		Qind.sim <- rnorm(n, mean=mQ, sd=sdQ)
		return(data.frame(Qind.sim=Qind.sim))
	}
	# Simulazione Qind regionale
	ARPIEM2012.sim.Qind_reg_moltiplicativo <- function(mQ, sdQ, n){
		CVQ <- sdQ/mQ
		teta1 <- log(mQ) - 0.5*log(1 + CVQ^2)
		teta2 <- sqrt(log(1 + CVQ^2))
		Qind.sim <- exp(rnorm(n, mean=teta1, sd=teta2))
		return(data.frame(Qind.sim=Qind.sim))
	}
	# Simulazione LCV e LCA entrambi campionari
	ARPIEM2012.sim.LCV_camp.LCA_camp <- function(mtau, sdtau, mtau3, sdtau3, n){
		tau3.sim <- rnorm(n, mean=mtau3, sd=sdtau3)
		ro <- (1-exp(-5*mtau3))/(1+exp(-5*mtau3))
		mtau.cond <- mtau + ro*sdtau * (tau3.sim - mtau3)/sdtau3
		sdtau.cond <- sdtau * sqrt(1 - ro^2)
		tau.sim <- rnorm(n, mean=mtau.cond, sd=sdtau.cond)
		return(data.frame(tau.sim = tau.sim, tau3.sim = tau3.sim))
	}
	# Simulazione LCV campionario e LCA regionale
	ARPIEM2012.sim.LCV_camp.LCA_reg_additivo <- function(mtau, sdtau, mtau3, sdtau3, n){
		tau.sim <- rnorm(n, mean=mtau, sd=sdtau)
		tau3.sim <- rnorm(n, mean=mtau3, sd=sdtau3)
		# Correzione per tau3.sim > 0.94 e tau.sim > 0.99
		tau3.sim[tau3.sim > 0.94] <- 0.94
		tau.sim[tau.sim > 0.99] <- 0.99
		return(data.frame(tau.sim = tau.sim, tau3.sim = tau3.sim))	
	}
	# Simulazione LCV e LCA entrambi regionali
	ARPIEM2012.sim.LCV_reg_moltiplicativo.LCA_reg_additivo <- function(mtau, sdtau, mtau3, sdtau3, n){
		CVtau <- sdtau/mtau
		teta1 <- log(mtau) - 0.5*log(1 + CVtau^2)
		teta2 <- sqrt(log(1 + CVtau^2))
		tau.sim <- exp(rnorm(n, mean=teta1, sd=teta2))
		tau3.sim <- rnorm(n, mean=mtau3, sd=sdtau3)
		# Correzione per tau3.sim > 0.94 e tau.sim > 0.99
		tau3.sim[tau3.sim > 0.94] <- 0.94
		tau.sim[tau.sim > 0.99] <- 0.99
		return(data.frame(tau.sim = tau.sim, tau3.sim = tau3.sim))	
	}

	
	# genera la piena indice random
	if (Qind.type == "R"){
		Qind.sim <- ARPIEM2012.sim.Qind_reg_moltiplicativo(Qind, sdQind, n)		
	} else if (Qind.type == "C") {
		Qind.sim <- ARPIEM2012.sim.Qind_camp(Qind, sdQind, n)		
	} else {
		stop("Errore nella definizione del modello per la piena indice")
	}

	# genera LCV e LCA random
	if (LCV_LCA.type == "C_C"){
		x.temp <- ARPIEM2012.sim.LCV_camp.LCA_camp(LCV, sdLCV, LCA, sdLCA, n)
	} else if (LCV_LCA.type == "C_R") {
		x.temp <- ARPIEM2012.sim.LCV_camp.LCA_reg_additivo(LCV, sdLCV, LCA, sdLCA, n)
	} else if (LCV_LCA.type == "R_R") {
		x.temp <- ARPIEM2012.sim.LCV_reg_moltiplicativo.LCA_reg_additivo(LCV, sdLCV, LCA, sdLCA, n)
	} else {
		stop("Errore nella definizione del modello per LCV-LCA")
	}
		
	Lmom.sim <- cbind(Qind.sim=Qind.sim[,1], LCV.sim=x.temp[, "tau.sim"], LCA.sim = x.temp[, "tau3.sim"])
	return(Lmom.sim)	
}
##############################################################################################



##############################################################################################	
ARPIEM2012.freq <- function(Qind, sdQind, LCV, sdLCV, LCA, sdLCA, Qind.type, LCV_LCA.type, n=10000, Tr=c(20,50,100,200,500), conf.bands=c(.1,.2,.8,.9)){
	# calcolo dei parametri della curva di frequenza (lognormale a 3 parametri)
	parln3 <- unlist(par.lognorm(Qind, LCV*Qind, LCA))
	# calcolo dei quantili
	Freq <- 1-1/Tr
	Qln3 <- invF.lognorm(Freq, parln3[1], parln3[2], parln3[3])
	# calcolo curva di crescita
	KTln3 <- Qln3 / Qind
	
	# generazione Lmomenti random
	Lmom.sim <- ARPIEM2012.sim.Lmoments(Qind.type, LCV_LCA.type, Qind, sdQind, LCV, sdLCV, LCA, sdLCA, n)
				
	# calcolo dei parametri della LN3 per i dati simulati
	parln3.sim <- apply(Lmom.sim, 1, function(x){tryCatch(unlist(par.lognorm(x["Qind.sim"], x["LCV.sim"]*x["Qind.sim"], x["LCA.sim"])), error=function(e) {out<-c(NaN,NaN,NaN); return(out)})})
	# calcolo dei quantili (ogni colonna è una curva simulata)
	qln3.sim <- apply(parln3.sim, 2, function(x){tryCatch(invF.lognorm(Freq, x[1], x[2], x[3]), error=function(e) {rep(NaN,length(Freq))})})
		
	# determino le fasce di confidenza
	qln3.bands <- apply(qln3.sim, 1, quantile, conf.bands, na.rm=TRUE)
	colnames(qln3.bands) <- paste("T.", Tr, sep="")
	
	output <- list(Freq = Freq,
					Tr = Tr,
					KTln3 = KTln3,
					qln3 = Qln3,
					qln3.bands = qln3.bands
					)
	return(output)
}
##############################################################################################	
