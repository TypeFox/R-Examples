solarangle <-
function(lat,Jday){
# angle of solar inclination from horizontal at solar noon [rad]

#lat: latitdue [rad]
#Jday: Julian date or day of the year [day]

# solar declination [rad]
dec<-declination(Jday)

return(asin(sin(lat)*sin(dec)+cos(lat)*cos(dec)*cos(0)))
}

