# LatticeKrig
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html


# test of radial basis function based on Wendland
# and using sparse formats
# Important check is of the FORTRAN function dfind2d
# that does pairwise distances among points within a specified range.

library(LatticeKrig)
options( echo=FALSE)
 test.for.zero.flag<- 1
set.seed(123)
x1<-  matrix( runif(3*5), ncol=3)
x2<-  matrix( runif(3*6), ncol=3)
n1<- nrow(x1)
n2<- nrow(x2)



# 3d
look1<-LKrigDistance(x1, x2, delta=.7)
look1<- spind2full(look1)
look2<- rdist( x1,x2)
look2[ look2>.7] <-0
test.for.zero( look1,look2)
# 2d
look1<-LKrigDistance(x1[,1:2], x2[,1:2], delta=.7)
look1<- spind2full(look1)
look2<- rdist( x1[,1:2],x2[,1:2])
look2[ look2>.7] <-0
test.for.zero( look1,look2)
################# more calls to distance method
test.for.zero.flag<- 1
set.seed( 123)
N<- 10
x1<- matrix( runif( 3*N), N,3)*4 + 1
gList<- list(1:3, 1:5, 1:4)
class( gList)<- "gridList"  

x2<- make.surface.grid( gList)
out1<- LKrigDistance(x1, x2, delta = 2 )
out2<- LKrigDistance(x1, gList, delta= 2)

test.for.zero( out1$ind, out2$ind, relative=FALSE,
 tag="agreement (ind) between calls with matrix and gridList")
test.for.zero( max(abs(out1$ra-out2$ra)),0, relative=FALSE,tol= 5e-7,
 tag="agreement (ra) between calls with matrix and gridList")
 out3<- rdist( x1, x2)
 out3[ out3 > 2] <- 0
test.for.zero( max(abs(spind2full(out1)-out3)),0, relative=FALSE,tol= 5e-7,
 tag="agreement (ra) between calls with rdist and gridList")

out4<- LKrigDistance(x1, gList, delta = 2, components=TRUE )
out5<- LKrigDistance( x1, x2, delta=2, components=TRUE)
#out6<- LKDistComponents( x1,x2, delta=2)
test.for.zero( out4$ind, out5$ind, relative=FALSE,
 tag="agreement (inf) between calls with matrix and gridList" )
test.for.zero( max(abs(out4$ra-out5$ra)),0, relative=FALSE,tol=5e-7,
 tag="agreement (ra) between calls with matrix and gridList" )

############## checking different metrics

# 
data(NorthAmericanRainfall)
x0<-  cbind( NorthAmericanRainfall$longitude, NorthAmericanRainfall$latitude)
x1<- x0[1:10,]  
x2<- x0[11:18,]

look2<- rdist( directionCosines(x1),directionCosines(x2))*4000
look2[ look2>100] <-0

dtype<- "Chordal"
attr(dtype, "Radius")<- 4000
look1<-LKrigDistance(x1, x2, delta= 100, distance.type=dtype)
look1<- spind2full(look1)
test.for.zero( look1,look2,tag="Chordal distance")

look1<-LKrigDistance(x1, x2, delta= 100, distance.type="GreatCircle")
look1<- spind2full(look1)
look2<- rdist.earth( x1,x2)
look2[ look2>100] <-0

test.for.zero( look1,look2,tag="Great Circle")

set.seed(123)
x1<-  matrix( runif(2*5), ncol=2)
x2<-  matrix( runif(2*6), ncol=2)
n1<- nrow(x1)
n2<- nrow(x2)

look1<-Radial.basis(x1,x2, basis.delta=.5)
look1<- as.matrix(look1)
look2<- Wendland2.2(rdist( x1,x2)/.5)
test.for.zero( look1,look2, tag="Radial.basis verses rdist")

#### check marginal variances this is an global test of basis function function
   data( ozone2)
   x<-ozone2$lon.lat
   y<- ozone2$y[16,]
   good <-  !is.na( y)
   x<- x[good,]
   y<- y[good]
   a.wght<- 5
   x<- x[1:10,]
   y<- y[1:10]
   LKinfo<- LKrigSetup(x, NC=4,  a.wght=a.wght, alpha=1, nlevel=1,
        normalize=FALSE, NC.buffer=1)
   xg<- make.surface.grid( list(x= seq( -90, -85,, 4), y= seq( 38, 42,,3)) )
   PHI1<- LKrig.basis(x, LKinfo)
   look1<- as.matrix( PHI1)
   dtemp<- LKinfo$latticeInfo$delta[1]* LKinfo$basisInfo$overlap
   # NOTE: 2d rectangle model returnd a gridList object for the 
   # centers so need to expand this into a grid of locations
   centerLocations<- 
         make.surface.grid( LKrigLatticeCenters(LKinfo, Level=1) )       
   look2 <-  WendlandFunction( rdist(x, centerLocations)/dtemp )
   test.for.zero( look1, look2, relative=FALSE, tol= 4e-6,
    tag="check basis functions with Euclidean distance compare to rdist")
   PHI1test<-  Radial.basis( x, LKrigLatticeCenters(LKinfo, Level=1),
        basis.delta = LKinfo$latticeInfo$delta[1]*LKinfo$basisInfo$overlap)
   test.for.zero( PHI1, PHI1test, tag="check basis functions Euclidean distance")     
   PHI2<- LKrig.basis(xg, LKinfo)                  
   Q<- LKrig.precision( LKinfo)
   Ktest1<- PHI1%*%solve(Q)%*%t(PHI2)
   test.for.zero( Ktest1, LKrig.cov( x,xg, LKinfo=LKinfo))
   Ktest2<- PHI2%*%solve(Q)%*%t(PHI2)
   test.for.zero( diag(Ktest2), LKrig.cov( xg, LKinfo=LKinfo, marginal =TRUE),
                   tag="marginal variance")      
#                               
# lon lat geometry but with Great Circle distance
#
   LKinfo<- LKrigSetup(x, NC=4,  a.wght=a.wght, alpha=1, nlevel=1,
                           normalize = FALSE,
                       distance.type = "GreatCircle",
                           NC.buffer = 1) 
# artifically reset lattice spacing to be on the order of GC distance                          
                           LKinfo$latticeInfo$delta<- 100      
           PHI1<- LKrig.basis(x, LKinfo)
           PHI1test<-  Radial.basis( x, LKrigLatticeCenters(LKinfo, Level=1),
        basis.delta = LKinfo$latticeInfo$delta[1]*LKinfo$basisInfo$overlap,
        distance.type = LKinfo$distance.type)
       test.for.zero( PHI1, PHI1test, tag="check basis functions great circle")
#
cat( "Done with testing LKrig basis", fill=TRUE)
options( echo=TRUE)
