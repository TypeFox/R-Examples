ElliView <-
function(elli.dat, out.file, obj.res = 3){
  #Determine angle of rotation about the X axis
  x.rot <- ifelse(elli.dat$x.az <= 180,
  ifelse(elli.dat$y.az > elli.dat$x.az,
    ifelse(elli.dat$y.az <= (elli.dat$x.az + 180),
    elli.dat$y.in / cos(elli.dat$x.in * (pi / 180)),
    elli.dat$y.in / cos(elli.dat$x.in * (pi / 180)) * (-1)),
  elli.dat$y.in / cos(elli.dat$x.in * (pi / 180)) * (-1)),
  ifelse(elli.dat$y.az < elli.dat$x.az,
  ifelse(elli.dat$y.az <= (elli.dat$x.az + 180),
    elli.dat$y.in / cos(elli.dat$x.in * (pi / 180)) * (-1),
    elli.dat$y.in / cos(elli.dat$x.in * (pi / 180))),
  elli.dat$y.in / cos(elli.dat$x.in * (pi / 180)))
  )
  
  #Determine Y and Z axis rotation
  y.rot <- elli.dat$x.in
  z.rot <- elli.dat$x.az + 90
  
  #Create sphereical triangular mesh to distort into ellipsoid
  my.sphere <- subdivision3d(icosahedron3d(col = "gray"), depth = obj.res)
  
  #Distort and rotate sphere into fabric ellipsoid
  x.rot <- rotationMatrix((x.rot * (pi / 180)), 1, 0, 0)
  y.rot <- rotationMatrix((y.rot * (pi / 180)), 0, 1, 0)
  z.rot <- rotationMatrix((z.rot * (pi / 180)), 0, 0, 1)
  my.ellip <- scale3d(my.sphere, elli.dat$x.mag, elli.dat$y.mag, elli.dat$z.mag)
  my.ellip <- rotate3d(obj = my.ellip, matrix = x.rot)
  my.ellip <- rotate3d(obj = my.ellip, matrix = y.rot)
  my.ellip <- rotate3d(obj = my.ellip, matrix = z.rot)
  
  #Draw ellipsoid scene
  open3d()
  shade3d(my.ellip)
  
  #Write object file if desired
  if(!missing(out.file)){
    setwd("output")
    writePLY(out.file)
setwd("..")
  }
  
  #Add reference frame to scene
  box3d(col="#DDDDDD",lwd=.5)
  mtext3d("N", edge='x++',col="#DDDDDD")
  mtext3d("S", edge='x-+',col="#DDDDDD")
  mtext3d("E", edge='y++',col="#DDDDDD")
  mtext3d("W", edge='y-+',col="#DDDDDD")
}
