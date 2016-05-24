reflectance <- function(date,latitude,longitude,timezone){
  
  #Calculate sunvector for each date/time
  sv <- sunvector(JD(date),latitude,longitude,timezone)
  
  #Calculate corresponding zenith angle
  zenith <- sunpos(sv)[,2]
  
  #Reflectance is calculated as a function of zenith angle following Kirk 2011
  zen <- c(seq(0,85,5),87.5,89) 
  r   <- c(0.020,0.020,0.020,0.020,0.020,0.021,0.021,0.022,0.024,0.028,
           0.033,0.043,0.059,0.086,0.133,0.211,0.347,0.583,0.761,0.896)
  
  reflectance <- approx(zen,r,xout=zenith,rule=2)$y
  
  
  #Compute decimal day
  decday <- strptime(date,format="%Y-%m-%d %H:%M")$yday + 
            strptime(date,format="%Y-%m-%d %H:%M")$hour/24 + 
            strptime(date,format="%Y-%m-%d %H:%M")$min/24/60
  
  return(cbind(decday, reflectance))
  
}
