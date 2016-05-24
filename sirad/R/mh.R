mh <-
  function(days,lat,Tmax,Tmin) {
    i <- dayOfYear(days)
    latt <- radians(lat)
    LD <- 180/pi*0.267*asin((0.5+0.007895/cos(latt)+0.2168875*tan(latt))^0.5)   # longest day
    MHa <- (sin(latt)*(46.355*LD-574.3885)+816.41*cos(latt)*sin(pi*LD/24))*(0.29*cos(latt)+0.52)
    MHb <- (sin(latt)*(574.355-1.509*LD)-26.59*cos(latt)*sin(pi*LD/24))*(0.29*cos(latt)+0.52)
    Is <- 0.04188*(MHa+MHb*sin(((2*pi*(i+10.5))/365)-(pi/2)))
    tal <- 0.8+0.12*(abs(182-i)/183)^1.5
    mh <- (0.182*(Tmax-Tmin)^0.69*(tal*Is)^0.91-2.4999)/0.8023
    mh 
  }

