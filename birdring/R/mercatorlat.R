# computes the factor of the latitudes for the mercator projection
# Rummler 2002
# Author: Fraenzi Korner, 2004, www.oikostat.ch
mercatorlat<-function(x) {
##Funktion aus Rummler 2002, Elem. Math.
log(tan(pi/4+x*pi/360))
}
