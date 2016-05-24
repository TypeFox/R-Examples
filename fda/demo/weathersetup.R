#  Sets up the daily weather data as a named list containing:
#  weather$temp:  a 365 by 35 matrix containing average temperature 
#         in degrees Celsius for
#  weather$prec:  a 365 by 35 matrix containing average precipitation 
#         in millimetres for
#  weather$time:  a vector of 365 mid-day times in days
#  weather$names: a vector of station names

#  The list is saved with file name "weatherdata" and
#  can be loaded with the command
#  load("weatherdata")

#  These commands are executed in folder fdaR/demo

#  Last modified 17 November 200

#  ------------------------  set up the data  -----------------------

tempav  <- daily$tempav
precav  <- daily$precav
station <- daily$place

#  set up the times of observation at noon

daytime   <- (1:365)-0.5

weatherdata <- list(tempav  = tempav,  precav = precav,
                    daytime = daytime, station = station)

save(weatherdata, file = "weatherdata")

