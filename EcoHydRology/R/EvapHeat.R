EvapHeat <- function (surftemp, airtemp, relativehumidity=NULL, Tn=NULL, wind=2) {
    
	##	surftemp:  			Temperature of surface 	[degrees C]
	##	airtemp:  			Temperature of air	 	[degrees C]	
	##	relativehumidity:  	between 0 - 1			[-]
	##	Tn 					minimum dailiy air temperature, assumed to be the dewpoint temperature [C] 
	##	surftemp:  			Temperature of surface 	[degrees C]
	## 	wind 				average daily windspeed [m/s] 
	
	
	windfunction <- 5.3 * (1 + wind)
    if (relativehumidity >= 0 & relativehumidity <= 1) {
        airvapordensity <- relativehumidity * SatVaporDensity(airtemp)
    }
    else {
        airvapordensity <- SatVaporDensity(Tn)
    }
    surfacevapordensity <- SatVaporDensity(surftemp)
    return(round(86400 * windfunction * (surfacevapordensity - airvapordensity)))
}

