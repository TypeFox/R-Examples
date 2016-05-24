
incident <- function(date, latitude, longitude, elevation,
                     timezone, meanPAR, TL=3.5, reflectance=TRUE){
  
#Calculate sunvector and sunposition for each date/time
sv <- sunvector(JD(date),latitude,longitude,timezone)

sp <- sunpos(sv)

######################
#Implement Hofierka and Suri 2002
######################

#Calculate Decimal Day of Year
decday <- strptime(date,format="%Y-%m-%d %H:%M")$yday + 
          strptime(date,format="%Y-%m-%d %H:%M")$hour/24 + 
          strptime(date,format="%Y-%m-%d %H:%M")$min/24/60

#Calculate solar constant (SC)
SC        <- 1367.13 * (1 + 0.03344 *cos(2*pi*decday/365.25-0.048869))

#Calculate solar altitude (ho) in degrees
ho <- 90-sunpos(sv)[,2]

#Calculate corrected solar altitude (horef) due to atomspheric refraction. 
#Eq. 6 in Hofierka and Suri
horef <- ho + 0.061359*(0.1594+1.123*ho+0.065656*ho*ho)/(1+28.9344*ho+277.3971*ho*ho)

#Negative solar altitudes are assigned 0
ho[ho<0]       <- 0
horef[horef<0] <- 0

#Calculate relative optical air mass (m). Eq 5 in Hofierka and Suri
m <- exp(-elevation/8434.5)/(sin(horef*pi/180)+0.50572*(horef+6.07995)^(-1.6364)) 

#Calculate Rayliegh optical thickness(delr). Eq 8-9 in Hofierka and Suri 
delr        <- 1/(6.6296 + 1.7513*m - 0.1202*m^2 + 0.0065*m^3- 0.00013*m^4)
delr[m>20]  <- 1/(10.4 + 0.718*m[m>20])

#Calculate beam irradiance (B)
B   <- SC * exp(-0.8662 * TL * m * delr)
B   <- B * sin(ho*pi/180)

#Calculate Diffuse radiation (D)
Tn <- -0.015843 + 0.0305843*TL - 0.0003797*TL*TL

#Coefficients. Eq 24 in Hofierka and Suri
A1 <- 0.26463 - 0.061581*TL + 0.0031408*TL*TL
A1[A1<0.0022]  <- 0.0022/Tn[A1<0.0022]
A2 <- 2.04020 + 0.018945*TL - 0.011161*TL*TL
A3 <- -1.3025 + 0.039231*TL + 0.0085079*TL*TL

#Calculate diffuse solar radiation function (Fd). Eq. 23 in Hofierka and Suri
Fd <- A1+A2*sin(ho*pi/180)+A3*sin(ho*pi/180)^2

D  <- SC*Tn*Fd
#Set nocturnal values to 0
D[ho==0] <- 0

#Total incident irradiance sum of beam and diffuse radiation
E0    <- B + D
#Convert shortwave radiation to PAR, W/m2 to umol m-2 s-1
E0    <-  E0 * 4.6 * 0.445  

######################
#Scale PAR to user data
######################

#Scale E0 to user value if it exists
if (missing("meanPAR")){

  #No value passed so do nothing
  
} else {  
  
  #Determine number of simulation days
  n  <- difftime(max(date),min(date),"days")
  n  <- as.numeric(n)
  
  #Integrate E0
  fn    <- splinefun(JD(date),E0,method="periodic")    
  meanE <- integrate(fn,min(JD(date)),max(JD(date)))$value / n
  
  #Scale
  E0 <- E0/meanE*meanPAR
  
}

#Subtract surface reflectance if desired following Kirk 2011
if (reflectance==TRUE){
  ref <- reflectance(date,latitude,longitude,timezone)
  E0 <- E0 - E0*ref[,2] 
}

return(cbind(decday, E0))

}
