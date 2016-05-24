pointsphere <-
function(N=100,longlim=c(0,360),latlim=c(-90,90),rlim=c(0,1)){
long=runif(N,longlim[1],longlim[2])
lat=asin(runif(N,sin(latlim[1]*pi/180),sin(latlim[2]*pi/180)))*180/pi
radius=runif(N,rlim[1]^3,rlim[2]^3)^{1/3}
return=cbind(long=long,lat=lat,radius=radius)
}
