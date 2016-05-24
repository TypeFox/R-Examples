declination <-
function(Jday){
# solar declination [rad]
#
#Jday: Julian date or day of the year [day]
#
return(0.4102*sin(pi*(Jday-80)/180))
}

