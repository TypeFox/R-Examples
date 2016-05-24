SM <- function(depth, min_interval, DO, temp, K, day, sr="00:00:00", ss="23:45:00", start="00:00:00", end="23:45:00"){

#######################################################################
#Stephen A. Sefick
#This function was written following:
#Izagirre, O., M. Bermejo, J. Pozo, and A. Elosegi. 2007. RIVERMET: An Excel-based tool to calculate river metabolism from diel oxygen concentration curves. Environmental Modelling and Software, 22: 24-32.
#depth m
#DO mg/L
#temp degrees Celcius
#K 1/d

#Output is in Oxygen mg/L*d

#######################################################################
  
#Dissolved Oxygen Deficit
D <- Cs(temp) - DO

#correct K20 to stream temp
Ktemp <- Kt(K, temp)

#Ktemp 0-sunrise; sunset-0
Ktemp.sr <- window_chron(Ktemp, day, start, day, sr) 
Ktemp.ss <- window_chron(Ktemp, day, ss, day, end)

#Rearation flux
Re_flux <- Ktemp*D/min_interval

#DO 0-sunrise; sunset-0
DO.sr <- window_chron(DO, day, start, day, sr) 
DO.ss <- window_chron(DO, day, ss, day, end)

#Respiration 0-sunrise; sunset-0
R.night.sr <- (dC.dt(DO.sr)-Ktemp.sr*D)/min_interval
R.night.ss <- (dC.dt(DO.ss)-Ktemp.ss*D)/min_interval

#mean nighttime Respiration
Rnight <- mean(coredata(rbind(R.night.sr, R.night.ss)), na.rm=T)

#Temperature 0-sunrise; sunset-0
temp.sr <- window_chron(temp, day, start, day, sr) 
temp.ss <- window_chron(temp, day, ss, day, end)

#Average Nightime Temperature
avg_night_temp <- mean(coredata(rbind(temp.sr, temp.ss)), na.rm=T)

#Daytime Temperature
day_temp <- window_chron(temp, day, sr, day, ss)

#Daytime Respiration
Rday <- Rnight*(1.072^(day_temp-avg_night_temp))

#Daytime DO
DO.day <- window_chron(DO, day, sr, day, ss)

#D in the daytime
D.day <- window_chron(D, day, sr, day, ss)
#################

#Community Respiration 24 hour corrected for time
ER <- Rnight*(1.072^(temp-avg_night_temp))*min_interval

#Net Ecosystem Production
NEP <- dC.dt(DO)-Ktemp*(Cs(temp)-DO)

NEP.ER <- merge(NEP, ER)


#GPP during daytime
GPP <- sum(window_chron(NEP.ER[,1]-NEP.ER[,2], day, sr, day, ss))

#make into g/m_sq*d
GPP*depth
ER*depth
NEP*depth

#NEP is the sum of all day CR and light time GPP
return(data.frame(GPP, ER=sum(ER), NEP=sum(ER)+GPP))

}
