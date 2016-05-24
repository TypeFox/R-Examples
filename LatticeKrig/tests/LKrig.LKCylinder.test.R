# LatticeKrig
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library( LatticeKrig)
options( echo=FALSE)

##########################################
  test.for.zero.flag<- 1
    

################ check of basis functions using explicit formula
NC<- 4
xdomain<- cbind( c(0,360), c( -70,70), c( 1,2.7))
LKinfo<- LKrigSetup(
         xdomain,
           NC = NC,
       nlevel = 1,
   LKGeometry = "LKCylinder",
       a.wght = 6.2,
    NC.buffer = 0,
    normalize = FALSE,
    BasisFunction = "triWeight",
          overlap = 2.5,
          mean.neighbor= 400,
    choleskyMemory = list( nnzR = 5e7),
                 V = diag( c( 1, 1, (2.7-1)/180 ))
     )


  sphericalCenters<- LKrigLatticeCenters( LKinfo, Level=1, 
                          physicalCoordinates=TRUE)$Locations
  xyzCenters<- directionCosines( sphericalCenters )
#  library( rgl); plot3d(xyzCenters)



   sphericalCenters<- LKrigLatticeCenters( LKinfo, Level=1, 
                          physicalCoordinates=TRUE)$Locations
  NIndex<- 20
# look at basis function with index  NIndex
   nodeCenter<-  rbind( sphericalCenters[NIndex,] )
#
   a<- LKrigLatticeCenters( LKinfo, Level=1, 
                          physicalCoordinates=TRUE)$basisScaling 
# create a grid but not exactly at node points. 
   gridL<- list( seq(0,360,,45), seq( -60,60,,45),
                    nodeCenter[3] +.2  )
   xGrid<- make.surface.grid( gridL)
   N<- nrow( xGrid)                        
   value<- rep( NA, N) 
# handy distance function for angle in [0,360]    
   distance360<- function( u1,u2){
   	 	min( abs(u2-u1), abs(u1 - u2 + 360), abs( u2- u1 + 360) )
   	}
   
   for( k in 1:N){
# first coordinate is periodic [0,360]   	
     d1<- ( distance360( xGrid[k,1], nodeCenter[1] ) /a[1] )^2
   	 d2<- ( (xGrid[k,2]-nodeCenter[2]       )/a[2]  )^2
   	 d3<- ( (xGrid[k,3]-nodeCenter[3]       )/a[3]  )^2
     distance2<- (d1) + (d2)  + (d3)  
 # triweight    
   	 value[k]<- ifelse( distance2 < 1, (1 - distance2)^3, 0) 
   	}

# compare to more efficient computation within LatticeKrig    		
     bigX<- LKrig.basis( xGrid, LKinfo)
      
     value1<- bigX[,NIndex]

# maximum error    
    
     test.for.zero( value, value1, tol=1e-10,
      tag="LKCylinder basis using expliciti formula")  