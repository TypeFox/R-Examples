# examples for the mvmesh package

# 2d plots
plot( UnitSimplex( n=2, k=4 ), show.points=TRUE )
plot( SolidSimplex( n=2, k=4 ), col="cyan", show.faces=TRUE )
plot( UnitSphere( n=2, k=2 ), show.points=TRUE )
plot( UnitBall( n=2, k=2 ), show.faces=TRUE, col="green" )
plot( RectangularMesh( a=c(1,3), b=c(2,10), breaks=c(5,8) ) )
plot( PolarSphere( n=2), col="blue" )
plot( PolarBall( n=2 ), show.labels=TRUE )

# translate and rotate in 2d
plot(0,0,xlim=c(0,2),ylim=c(0,2),type='n') # set plot window
mesh <- SolidSimplex( n=2, k=3 )
plot(mesh, new.plot=FALSE, col="blue")
A <- Rotate2D( pi )
mesh2 <- AffineTransform( mesh, A=A, shift=c(1.1,1.1) ) 
plot(mesh2, new.plot=FALSE, col="red" )
mesh3 <- AffineTransform( mesh, A=diag( c(0.5,0.5) ), shift=c(1.2,1.2) )
plot(mesh3, new.plot=FALSE, col="green", show.faces=TRUE )


# different ways to specify subdivision of a rectangular mesh
plot( RectangularMesh( a=c(1,3), b=c(2,7) ) )
plot( RectangularMesh( a=c(1,3), b=c(2,7), breaks=c(4,10) ) )
plot( RectangularMesh( a=c(1,3), b=c(2,7), 
    breaks=list( seq(1,3,by=0.25), seq(2,7,by=1) ) ) )


#  3d plots
plot( UnitSimplex( n=3, k=4 ), col="magenta" )
plot( UnitSimplex( n=3, k=4 ), col="magenta", show.faces=TRUE )

plot( SolidSimplex( n=3, k=4 ), col="orange" )
plot( SolidSimplex( n=3, k=4 ), col="orange", show.faces=TRUE, alpha=0.5 )

plot( UnitSphere( n=3, k=2 ), col="blue")
plot( UnitSphere( n=3, k=2 ), col="blue", show.faces=TRUE )

plot( UnitBall( n=3, k=2 ), col="green" )
plot( UnitBall( n=3, k=2 ), col="green", show.faces=TRUE, alpha=0.5 )

plot( RectangularMesh( a=c(1,3,0), b=c(2,5,2), breaks=c(2,3,4) ), show.labels=TRUE, col="red" )
plot( RectangularMesh( a=c(1,3,0), b=c(2,5,2), breaks=c(2,3,4) ), show.labels=TRUE, col="red", show.faces=TRUE, alpha=0.5 )

plot( PolarSphere( n=3 ), col="blue")
plot( PolarSphere( n=3 ), col="blue",show.faces=TRUE )

plot( PolarBall( n=3 ), col="red" )

plot( Icosahedron( ), col="cyan" )
plot( Icosahedron( ), col="cyan", show.faces=TRUE )

# unit simplex with different subdivisions
A <- diag( rep(1.0,3) )
my.color <- c("black","red","blue","green")
for (k in 1:4) {
  mesh <- UnitSimplex( n=3, k=k )
  plot( AffineTransform( mesh, A=A, shift=c(k,0,0)), new.plot=(k==1),col=my.color[k] )
  text3d( k,0,1, paste("k=",k) )
}
title3d("different subdivions")


# translate and rotate in 3d
mesh <- SolidSimplex( n=3, k=2 )
plot(mesh, col="blue")
A <- Rotate3D( rep(pi/2,3) )
mesh2 <- AffineTransform( mesh, A=A, shift=c(1,1,1) ) 
plot(mesh2, new.plot=FALSE, col="red" )

# different types of tubes
p.range <- c(0.5,1,2,4,Inf)
open3d()
mfrow3d( nr=2, nc=length(p.range), byrow=FALSE, sharedMouse=TRUE )
for (i in 1:length(p.range)) {
  next3d()
  plot( HollowTube( n=3, k.x=3, p=p.range[i]), new.plot=FALSE, col=i ) 
  next3d()
  plot( SolidTube( n=3, k.x=3, p=p.range[i]), new.plot=FALSE, col=i )
  title3d( paste("tubes with p=",p.range[i]) )
}


# check UnitSimplex in different dimensions
for (n in 2:5) {
  for (k in 1:4) {
    mesh <- UnitSimplex( n=n, k=k )
    cat("UnitSimplex in dimension n=",n,"  k=",k,"  vertices=",ncol(mesh$V),"  simplices=",dim(mesh$S)[3],"\n")
  }
}



