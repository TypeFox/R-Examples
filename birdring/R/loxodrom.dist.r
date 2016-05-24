################################################################################
## Computes the distance along the loxodrome in km between two points on the Earth
## x1, x2: longitudes of points 1 and 2 in decimal coordinates, W = negative sign
## y1, y2: latitudes  of points 1 and 2 in decimal coordinates, S = negative sign
## epsilon: threshold value of a number to be interpreted as zero
## package: either "geosphere" or "birding". If "birdring" the "old" function written by F Korner is used. This function does not work in every case... the new function uses the function distRhomb from the package geosphere which is much better! 
## Reference: Imboden, C. & D. Imboden (1972). Vogelwarte 26: 336-346.
## Author: Fraenzi Korner, Sept. 2004, www.oikostat.ch/www.vogelwarte.ch
## depends on package geosphere
#################################################################################
loxodrom.dist<-function(x1, y1, x2, y2, epsilon=0.0001, package="geosphere"){
if(package=="geosphere"){
  dist <- distRhumb(cbind(x1, y1), cbind(x2, y2))/1000
  return(dist)
}# close geosphere
  
# old version  
if(package=="birdring"){
dis<-numeric(length(x1))
rerde<-6368
deltax<-abs(x2*pi/180-x1*pi/180)
deltay<-abs(y2*pi/180-y1*pi/180)
tga<-deltax/(log(tan(pi/4+y2*pi/360))-log(tan(pi/4+y1*pi/360))) 

dis[abs(x1-x2)<epsilon&abs(y1-y2)<epsilon]<-0
dis[abs(y1-y2)<epsilon&(abs(x1-x2)>epsilon)]<-abs(cos(y1[abs(y1-y2)<epsilon&(abs(x1-x2)>epsilon)]*pi/180)*deltax[abs(y1-y2)<epsilon&(abs(x1-x2)>epsilon)])
dis[(tga<0)&(abs(x1-x2)>epsilon)&(abs(y1-y2)>epsilon)]<-abs(deltay[(tga<0)&(abs(x1-x2)>epsilon)&(abs(y1-y2)>epsilon)]/cos((pi-atan(tga[(tga<0)&(abs(x1-x2)>epsilon)&(abs(y1-y2)>epsilon)]))))
dis[(tga>=0)&(abs(x1-x2)>epsilon)&(y1-y2>epsilon)]<--deltay[(tga>=0)&(abs(x1-x2)>epsilon)&(y1-y2>epsilon)]/cos(atan(tga[(tga>=0)&(abs(x1-x2)>epsilon)&(y1-y2>epsilon)]))
dis[(tga>=0)&(abs(x1-x2)>epsilon)&(y2-y1>epsilon)]<-abs(deltay[(tga>=0)&(abs(x1-x2)>epsilon)&(y2-y1>epsilon)]/cos(atan(tga[(tga>=0)&(abs(x1-x2)>epsilon)&(y2-y1>epsilon)])))
dis[(abs(x1-x2)<epsilon)&(y2-y1>epsilon)]<-abs(deltay[(abs(x1-x2)<epsilon)&(y2-y1>epsilon)]/cos(atan(tga[(abs(x1-x2)<epsilon)&(y2-y1>epsilon)])))
dis[(abs(x1-x2)<epsilon)&(y1-y2>epsilon)]<-abs(deltay[(abs(x1-x2)<epsilon)&(y1-y2>epsilon)]/cos(atan(tga[(abs(x1-x2)<epsilon)&(y1-y2>epsilon)])))
return(dis*rerde)
}# close package birdring
}
###################################################################################