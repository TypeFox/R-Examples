

#########################################################
## r
## based on Schoolfield, R. M., Sharpe, P. J. & Magnuson, C. E. 
## 1981. Non-linear regression of biological temperature-dependent 
## rate models based on absolute reaction-rate theory. Journal of 
## Theoretical Biology, 88, 719-731.
##########################################################

.SSM<-function (T, parms) {

  # Si je travaille en degree celsius, je convertis en Kelvin
  if (T[1]<273) T <- T+273.15
  
	if (all(names(parms)!="Rho25")) {
	  newx <- abs(parms[(names(parms)!="rK") & (names(parms)!="K") & (names(parms)!="Scale")])
    scale <- ifelse(is.na(parms["Scale"]), 1, parms["Scale"])
#	  tableT <- data.frame(Temperature=as.numeric(names(newx)), R=newx)
#	  ml <- loess(R ~ Temperature, data=tableT)
#	  newtable<- data.frame(Temperature=T)
#	  rT <- predict(ml, newtable)*1E-5
    if (as.numeric(names(newx[1]))<273) names(newx) <- as.character(273.15+as.numeric(names(newx)))

#  	nnewx <- as.numeric(names(newx))
#		p <- poly.calc(x=nnewx, y=newx)
# 		r <- predict(p, T)

#	  rT <- ifelse(r<0, 0, r)*scale*1E-5
#	  rT_L <- rT

    nnewx <- as.numeric(names(newx))
	listpolynom <- NULL
	for (i in 1:(length(newx)-3)) {
		listpolynom <- c(listpolynom, list(list(
		  polynom=polynom::poly.calc(nnewx[i:(i+3)], newx[i:(i+3)]), 
		  range=c(nnewx[i+1], nnewx[i+2]))))
	}
	
	listpolynom[[1]]$range[1] <- nnewx[1]-10
	listpolynom[[length(listpolynom)]]$range[2] <- nnewx[length(nnewx)]+10

	r <- NULL

	for (i in 1:length(T)) {
		l <- unlist(lapply(listpolynom, function(p) ((T[i]>=p$range[1]) & (T[i]<=p$range[2]))))
 		pl <- which(l)[1]
 		p <- listpolynom[[pl]]$polynom
 		r <- c(r, predict(p, T[i]))
	}


	  rT <- ifelse(r<0, 1E-4, r)*scale*1E-5
	  rT_L <- rT
	
	} else {

    
	R <- 8.314472
	rho25 <- parms["Rho25"]/1E7
	dha <- parms["DHA"]*1E3
	dhh <- parms["DHH"]*1E3
  unsurT <- 1/T
  unsur298 <- 1/298

	if (is.na(parms["DHL"])) {
# Je suis en 4 parametres
#	"T12H", "DHA",  "DHH", "Rho25"
		t12H=abs(parms["T12H"])
		rT<-(rho25*(T/298)*exp((dha/R)*(unsur298-unsurT)))/(1+exp((dhh/R)*((1/t12H)-unsurT)))
	} else {
# 28/7/2012 - T12H change en DT	
#	"T12L", "DT", "DHA",  "DHH", "DHL", "Rho25"
		dhl <- parms["DHL"]*1E3
		t12L <- parms["T12L"]
		t12H <- t12L+abs(parms["DT"])
		rT<-(rho25*(T/298)*exp((dha/R)*(unsur298-unsurT)))/(1+exp((dhl/R)*((1/t12L)-unsurT))+exp((dhh/R)*((1/t12H)-unsurT)))
	}

	if (!is.na(parms["Rho25_L"])) {
	
	rho25_L <- parms["Rho25_L"]/1E7
	dha_L <- parms["DHA_L"]*1E3
	dhh_L <- parms["DHH_L"]*1E3

	if (is.na(parms["DHL_L"])) {
#	"T12H", "DHA",  "DHH", "Rho25"
		t12H_L <- abs(parms["T12H_L"])
		rT_L <- (rho25_L*(T/298)*exp((dha_L/R)*(unsur298-unsurT)))/(1+exp((dhh_L/R)*((1/t12H_L)-unsurT)))
	} else {
# 28/7/2012 - T12H change en DT	
#	"T12L", "T12H", "DHA",  "DHH", "DHL", "Rho25"
		dhl_L <- parms["DHL_L"]*1E3
		t12L_L <- parms["T12L_L"]
		t12H_L <- t12L+abs(parms["DT_L"])
		rT_L <- (rho25_L*(T/298)*exp((dha_L/R)*(unsur298-unsurT)))/(1+exp((dhl_L/R)*((1/t12L_L)-unsurT))+exp((dhh_L/R)*((1/t12H_L)-unsurT)))
	}
	} else {
	  rT_L <- rT
  }
	}

return(list(as.numeric(rT), as.numeric(rT_L)))

}
