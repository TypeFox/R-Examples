# show different types of meshes
UnitSimplex( n=3, k=4 )
SolidSimplex( n=3, k=4 )
UnitSphere( n=3, k=4, method="edgewise" )
UnitSphere( n=3, k=2, method="dyadic" )
UnitBall( n=3, k=3, method="edgewise" )
UnitBall( n=3, k=2, method="dyadic" )
RectangularMesh( a=c(1,3), b=c(2,10) )

# edge subdivision
T <- EdgeSubdivision( n=2, k=2 )
T
ConvertBase( 10, 2, 6 )  # note order of digits
NumVertices( n=4, k=8, single=FALSE )

S <- rbind( diag(rep(1,2)), c(0,0) ) # solid simplex in 2D
PointCoord( S, T[,,1] )
SimplexCoord( S, T[,,1] )
SVIFromColor( S, T )

# check all meshes in different dimensions
for (n in 2:3) {
  cat("dimension n=",n," ... ")
  timer <- system.time({
  UnitSimplex( n=n, k=4 )
  SolidSimplex( n=n, k=4 )
  UnitSphere( n=n, k=4, method="edgewise" )
  UnitSphere( n=n, k=2, method="dyadic" )
  UnitBall( n=n, k=3, method="edgewise" )
  UnitBall( n=n, k=2, method="dyadic" )
  RectangularMesh( a=rep(0.0,n), b=rep(1.0,n) ) })
  cat(timer,"\n")
}
