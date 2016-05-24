EnvirEnergy <-
function(lat,Jday,Tx,Tn,wind,rain,relativehumidity,albedo=0.2,cloudiness=NULL,forest=0,slope=0,aspect=0,surftemp=(Tx+Tn)/2,surfemissivity=0.97){
# the total energy exchange between the surface and the surrounding air

#lat: latitdue [rad]
#Jday: Julian date or day of the year [day]
#Tx: maximum daily temperature [C]
#Tn: minimum daily temperature [C]
#wind: average daily windspeed [m/s]
#relativehumidity: relative humidity; if negative, air vapor density will be approximated [-]
#cloudiness: fraction of the sky covered in clouds. If not provided, cloudiness will be approximated [-]
#albedo: surface albedo or reflectivity [-]
#forest: forest or vegeation cover [-]
#slope: slope of the ground [rad]
#aspect: ground aspect [rad from north]
#surftemp: surface temperature [C]
#surfemissivity: [-]
#rain: precipitation [mm/day]

if(is.null(cloudiness)){cloudiness<-EstCloudiness(Tx,Tn)}

airtemp<-(Tx+Tn)/2 #average daily air temperature [C]

return(Solar(lat,Jday,Tx,Tn,albedo,forest,slope,aspect)+Longwave(AtmosphericEmissivity(airtemp,cloudiness),airtemp)-Longwave(surfemissivity,surftemp)+SensibleHeat(surftemp,airtemp,wind)+EvapHeat(surftemp,airtemp,relativehumidity,Tn,wind)+RainHeat(airtemp,rain)+GroundHeat())
}

