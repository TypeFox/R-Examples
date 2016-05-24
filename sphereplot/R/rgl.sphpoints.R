rgl.sphpoints <-
function(long,lat,radius,deg=TRUE,col='black',...){
points3d(sph2car(long,lat,radius,deg=deg),col=col,...)
}
