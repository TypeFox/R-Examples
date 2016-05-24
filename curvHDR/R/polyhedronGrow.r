########## R function: polyhedronGrow ##########

# For growing a convex polyhedron by placing
# a ball of given radius "r" (in standard
# deviation units) in the center of each
# triangular face.

# Last changed: 03 APR 2009

polyhedronGrow <- function(cvxPolyhedron,r)
{
   Vertices <- unique(rbind(cvxPolyhedron$v1,cvxPolyhedron$v2,cvxPolyhedron$v3))

   stDevs <- apply(Vertices,2,sd)

   if (any(stDevs==0))
      stop("zero standard deviations inside polyhedronGrow().")
   
   sVertices <- Vertices%*%diag(1/stDevs)

   sConvexHull <- convhulln(sVertices,options="QJ")
   
   sv1 <- sVertices[sConvexHull[,1],]
   sv2 <- sVertices[sConvexHull[,2],]
   sv3 <- sVertices[sConvexHull[,3],]

   scvxPolyhedron <- makeTriangles(sv1,sv2,sv3,color="Red")
 
   currCentroid <- triPolyhCentroid(scvxPolyhedron)
 
   numFaces <- nrow(scvxPolyhedron$v1)

   sGrownPolyhPoints <- matrix(NA,numFaces,3)

   for (i in 1:numFaces)
   {
      # 3 vertices of the triangle:

      v1 <- scvxPolyhedron$v1[i,]
      v2 <- scvxPolyhedron$v2[i,]
      v3 <- scvxPolyhedron$v3[i,]

      # Set the rotation angles:
             			
      theta1 <- atan((v2[2]-v1[2])/(v2[1]-v1[1]))
      if ((v2[2]==v1[2])&(v2[1]==v1[1]))
         theta1 <- 0

      theta2 <- atan((v2[3]-v1[3])/(cos(theta1)*(v2[1]-v1[1])
                                    +sin(theta1)*(v2[2]-v1[2])))

      if ((v2[3]==v1[3])&(cos(theta1)*(v2[1]-v1[1])+sin(theta1)*(v2[2]-v1[2])==0))
         theta2 <- 0
 
      theta3 <- atan((-sin(theta2)*(cos(theta1)*(v3[1]-v1[1])
                     +sin(theta1)*(v3[2]-v1[2]))+cos(theta2)*(v3[3]-v1[3]))/
                     (-sin(theta1)*(v3[1]-v1[1])+cos(theta1)*(v3[2]-v1[2])))

      if ((-sin(theta2)*(cos(theta1)*(v3[1]-v1[1])+sin(theta1)*(v3[2]-v1[2]))
             +cos(theta2)*(v3[3]-v1[3])==0)
             &(-sin(theta1)*(v3[1]-v1[1])+cos(theta1)*(v3[2]-v1[2])==0))
         theta3 <- 0
         
      # Transform candidate points back to the original coordinates:
         
      q1back <- c((1/3)*(v1[1]+v2[1]+v3[1])+2*r*(sin(theta3)*sin(theta1)
                -cos(theta3)*sin(theta2)*cos(theta1)),(1/3)*(v1[2]+v2[2]+v3[2])
                 +2*r*(-sin(theta3)*cos(theta1)-cos(theta3)*sin(theta2)*sin(theta1)),
                 (1/3)*(v1[3]+v2[3]+v3[3])+2*r*cos(theta3)*cos(theta2))
 
      q2back <- c((1/3)*(v1[1]+v2[1]+v3[1])-2*r*(sin(theta3)*sin(theta1)
                 -cos(theta3)*sin(theta2)*cos(theta1)),(1/3)*(v1[2]+v2[2]+v3[2])
                 -2*r*(-sin(theta3)*cos(theta1)-cos(theta3)*sin(theta2)*sin(theta1)),
                 (1/3)*(v1[3]+v2[3]+v3[3])-2*r*cos(theta3)*cos(theta2))

      distPlus <- sqrt(sum((q1back-currCentroid)^2))
      distMnus <- sqrt(sum((q2back-currCentroid)^2))

      qback <- q1back
      if (distMnus>distPlus) qback <- q2back
         
      sGrownPolyhPoints[i,] <- qback
   }
 
   # Transform back to original scale:

   GrownPolyhPoints <- sGrownPolyhPoints%*%diag(stDevs)
 
   # Form the grown triangular polyhedron:

   gConvexHull <- convhulln(GrownPolyhPoints,options="FA")

   gv1 <- GrownPolyhPoints[gConvexHull$hull[,1],]
   gv2 <- GrownPolyhPoints[gConvexHull$hull[,2],]
   gv3 <- GrownPolyhPoints[gConvexHull$hull[,3],]

   grownCvxPolyhedron <- makeTriangles(gv1,gv2,gv3,color="green3",alpha=0.3)

   return(list(grownPolyHedron=grownCvxPolyhedron,volume=gConvexHull$vol))
}

########## End of polyhedronGrow ##########




