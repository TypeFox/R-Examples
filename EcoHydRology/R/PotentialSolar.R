PotentialSolar <-
function(lat,Jday){
# potential solar radiation at the edge of the atmospher [kJ m-2 d-1]

#lat: latitdue [rad]
#Jday: Julian date or day of the year [day]

# solar declination [rad]
dec<-declination(Jday)

return(117500*(acos(-tan(dec)*tan(lat))*sin(lat)*sin(dec)+cos(lat)*cos(dec)*sin(acos(tan(dec)*tan(lat))))/pi)
}

