plot3dInterval<-function(data,colors){
d1<-as.matrix(data[,,1])
d2<-as.matrix(data[,,2])
open3d(mouseMode="trackball")
#rgl.bg(sphere=T, texture=system.file("textures/sunsleep.png", package="rgl"), back="filled" )
for(i in 1:dim(data)[1]){
  #print(d1[i,])
  #print(d2[i,])
  t<-cube3d()
  wire3d(translate3d(scale3d(t,abs(d2[i,1]-d1[i,1])/2,abs(d2[i,2]-d1[i,2])/2,abs(d2[i,3]-d1[i,3])/2),(d2[i,1]+d1[i,1])/2,(d2[i,2]+d1[i,2])/2,(d2[i,3]+d1[i,3])/2),col=colors[i])
  #shade3d(scale3d(t,abs(d2[i,1]-d1[i,1]),abs(d2[i,2]-d1[i,2]),abs(d2[i,3]-d1[i,3])),col=col[i])
}
axes3d(edges = "bbox", labels = TRUE, tick = TRUE, nticks = 5)
#light3d(theta = 0, phi = 12)
TRUE
}

