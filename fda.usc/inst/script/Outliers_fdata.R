################################################################################
# Functional outlier detection procedure for detecting outliers via
# smooothed bootstrap procedure:
#  -based on trimming outliers.depth.trim()
#  -based on ponderation outliers.depth.pond()
#
# Three examples are shown in this script code using the datasets:
#  1. poblenou dataset
#  2. Canadian Weather dataset
#  3. aemet dataset
#
# For more details see, Febrero-Bande, M., Galeano, P. and Gonzalez-Manteiga, W.
# (2008).Outlier detection in functional data by depth measures with application
# to identify abnormal NOx levels}. Environmetrics 19, 4, 331{-}345.
#
################################################################################


################################################################################
#                        1.  poblenou dataset
#
# The aim is outlier detect  in the Nox levels split into two groups:
#(1) labour days.
#(2) nonworking days, which are the Saturdays, Sundays, and festive days.
################################################################################

data(poblenou)
names(poblenou)
# Functional data used
names(poblenou$nox)
# NOx levels measured every hour by a control station in Poblenou in Barcelona
# Each curve represents a day
nox<-poblenou$nox

# Non Functional data: poblenou$df
names(poblenou$df)
# $date
# $day.week: Factor levels:"Monday" 1, "Tuesday" 2, "Wednesday" 3, "Thursday" 4,
# "Friday" 5, "Saturday" 6 and "Sunday" 7.
# $day.festive: Factor levels: "non festive day" 0 and "festive day" 1.

working=poblenou$nox[poblenou$df$day.festive==0& as.integer(poblenou$df$day.week)<6]
nonworking=poblenou$nox[poblenou$df$day.festive==1 | as.integer(poblenou$df$day.week)>5]
out=outliers.depth.trim(nonworking,dfunc=depth.RP,nb=100,smo=0.1,trim=0.06)
out
out2= outliers.depth.trim(working,dfunc=depth.FM,nb=100,smo=0.1,trim=0.06)
out2

#Figure 3, Febrero et al (2008)
par(mfrow=c(2,1))
plot(nonworking,ylim=c(0,400),col="gray",lty=1)
lines(nonworking[out[[1]]],col=2,lty=2,lwd=2)
plot(working,ylim=c(0,400),col="gray",lty=1)
lines(working[out2[[1]]],col=2,lty=2,lwd=2)


################################################################################
#                        2.  Canadian Weather dataset
#
# The aim is outlier detect in the dayly temperature station
################################################################################
fdat<-fdata(t(CanadianWeather$dailyAv[,,1]))
m<-ncol(fdat)
n<-nrow(fdat)

# arguments
nb<-100;smo=0.05;trim=0.05
# It takes a lot (in funcion of number of resamples: nb )
#Method based on trimming: outliers.depth.trim
out.mode<-outliers.depth.trim(fdat,dfunc=depth.mode,nb=nb,smo=smo,trim=trim)
out.FM<-outliers.depth.trim(fdat,dfunc=depth.FM,nb=nb,smo=smo,trim=trim)
out.RP<-outliers.depth.trim(fdat,dfunc=depth.RP,nb=nb,smo=smo,trim=trim)
out.RPD<-outliers.depth.trim(fdat,dfunc=depth.RPD,nb=nb,smo=smo,trim=trim)

#Method based on ponderation: outliers.depth.pond
out2.mode<-outliers.depth.pond(fdat,dfunc=depth.mode,nb=nb,smo=smo)
out2.FM<-outliers.depth.pond(fdat,dfunc=depth.FM,nb=nb,smo=smo)
out2.RP<-outliers.depth.pond(fdat,dfunc=depth.RP,nb=nb,smo=smo)
out2.RPD<-outliers.depth.pond(fdat,dfunc=depth.RPD,nb=nb,smo=smo)

# plot
out<-out.FM
plot(fdat,col="gray",lty=1)
lines(fdat[out[[1]]],col=2,lty=2,lwd=2)



################################################################################
#                        3.  aemet dataset
#
# The aim is outlier detect  in the dayly temperature station
################################################################################
# Series of daily summaries of 73 spanish weather stations
data(aemet)
names(aemet)

# temperature curves by station
fdat<-aemet$temp
# arguments
nb<-10;smo=0.05;trim=0.05
# It takes a lot (in funcion of number of resamples: nb )
temp.mode<-outliers.depth.trim(fdat,dfunc=depth.mode,nb=nb,smo=smo,trim=trim)

# wind speed  curves by station
fdat<-aemet$wind.speed
# arguments
nb<-10;smo=0.1;trim=0.15
wind.FM<-outliers.depth.trim(fdat,dfunc=depth.FM,nb=nb,smo=smo,trim=trim)

# Log-precipitaton curves by station
fdat<-aemet$logprec
# arguments
nb<-10;smo=0.01;trim=0.15
prec.RP<-outliers.depth.trim(fdat,dfunc=depth.RP,nb=nb,smo=smo,trim=trim)

################################################################################
# Displaying the results

# Select an object to display the results
out<-temp.mode  # temperature curves
#out<-wind.FM    # wind speed curves
#out<-prec.RP   # log-precipitation curves

library(mapproj)
dev.new()
par(mfrow=c(1,3))

cl<-rainbow(length(out[[1]]))
#cl<-2:(length(out[[1]])-1)

# plot the curves
plot(fdat,col="gray",lty=1)
# plot the outliers detected
lines(fdat[out[[1]]],col=cl,lty=2,lwd=2)

# Representation on the map of the stations.
mapproject(map("world",c("Spain")),"mercator")
points(aemet$df$longitude,aemet$df$latitude,pch=1)
# Station with outliers in "Spain"
points(aemet$df$longitude[out[[1]]],aemet$df$latitude[out[[1]]],col=cl,pch=4,lwd=3)


mapproject(map("world",c("Canary Island")),"mercator")
points(aemet$df$longitude,aemet$df$latitude,pch=1)
# Station with outliers in "Canary Island"
points(aemet$df$longitude[out[[1]]],aemet$df$latitude[out[[1]]],col=cl,pch=4,lwd=3)



#mapproject(map("world",c("Spain","Canary Island")),"mercator")
#points(aemet$df$longitude,aemet$df$latitude,pch=1)
#points(aemet$df$longitude[out[[1]]],aemet$df$latitude[out[[1]]],col=cl,pch=4,lwd=3)



