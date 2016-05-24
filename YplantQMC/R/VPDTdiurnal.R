    

    VPDTdiurnal <- function(Tmax, Tmin,reltime=seq(0,1,length=12),
                        VPDmax=NA, percsatTmin=0.9, maxlag=0.1){

        # Saturated vapor pressure based on temperature.
        esatfun <- function(temp) 613.75*exp((17.502*temp)/(240.97 + temp))

    		# saturation at min. temperature.
    		eTmin <- percsatTmin * esatfun(Tmin)
    		if(is.na(VPDmax)){
    			VPDmax <- (esatfun(Tmax) - eTmin) / 1000
    		}
		
        # Temperature diurnal function
        # timerel : (0-1) relative time (during daylight hours)
        Tdiurnal <- function(timerel,Tmax,Tmin,maxlag){
            (Tmax-Tmin)*sin(pi*timerel / (1 + 2*maxlag)) + Tmin
        }

        # Vector of temperatures, from the diurnal function, evaluated at "dist".
        temperatures <- Tdiurnal(timerel=reltime,Tmax=Tmax,Tmin=Tmin,maxlag=maxlag)
					 
        # Vector of esat at these temperatures.
        esats <- esatfun(temperatures)

        # Precalculate VPD, then scale it to max and min VPD. 
        VPDs <- (esats - eTmin)/1000
		
        # Make sure not to mess with the minimum VPD. 
		    vpd1 <- VPDs - min(VPDs)
		    vpd2 <- vpd1 * (VPDmax-min(VPDs))/max(vpd1)
        VPDs <- vpd2 + min(VPDs)
		
    return(data.frame(reltime=reltime,VPD=VPDs, Tair=temperatures, RH=VPDtoRH(VPDs,temperatures)))
    }


