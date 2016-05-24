PotSolarInst <-
function(Jday, hour = 12, lat = 42.44*pi/180, sunrise = NULL, sunset = NULL, SolarNoon = mean(c(sunrise,sunset)), units = "Wm2", latUnits = "unknown"){
#  lat[assumes rad,  but can handle degrees if abs. value above 1.5]
#Jday[day of year]
#hour[hour of day, 0-24]
#Either sunrise and sunset times, or SolarNoon is needed   [in hours, 0-24]

if ((abs(lat) > pi/2 & latUnits == "unknown") | latUnits == "degrees" ){
lat <- lat*pi/180
} else if (latUnits == "unknown"){
warning("in PotSolarInst call: Latitude assumed to be in radians, if using degrees, please set latUnits = 'degrees'")
} 


SolarConstant <- 118000 #  [kJ/m2/d]
AngVeloc <- 0.2618  #rad/hr = 15 deg/hr
DayAngle <- 2*pi*(Jday-1)/365
dec <- declination(Jday)
ZenithAngle <- acos(sin(lat)*sin(dec) + cos(lat)*cos(dec)*cos(AngVeloc*(hour-SolarNoon)))
if (units == "Wm2") convert <- 86.4 else convert <- 1
PotSol <- SolarConstant * cos(ZenithAngle) / convert
PotSol[which(PotSol < 0)] <- 0
return ( signif(PotSol, 3) )
}
