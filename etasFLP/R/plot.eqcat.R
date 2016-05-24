plot.eqcat <-
function(x,...){
# defaul plots for earthquakes catalogs
# plot x-y points with magnitude extension and time colouring (red: recent, blu:older)

ts=(x$time-min(x$time))/diff(range(x$time))

mapxy=map("worldHires",xlim=range(x$long),ylim=range(x$lat))

points(x$long,x$lat,cex=sqrt(exp(x$magn1))/8,col=rgb(ts,0,1-ts),pch=19)
title(main="Circles extension proportional to magnitude \n red: recent, blu:older")


typegraph=2
plot3d(x$long,x$lat,x$time,type="n",zlab="time ",
xlab="x-longitude",ylab="y-latitude")
lines3d(cbind(mapxy$x,mapxy$y,min(x$time)),col="red")
lines3d(cbind(mapxy$x,mapxy$y,max(x$time)),col="red")
plot3d(x$long,x$lat,x$time,add=TRUE, type="p",col=rgb(ts,0,1-ts))

#####################################################
dev.new()
magn.plot(x)

}
