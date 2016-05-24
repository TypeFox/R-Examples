########## R function: ptsInTriPolyh ##########

# Computes indicators of whether points 
# are inside a convex triangular-faced polyhedron.

# Last changed: 17 JUL 2009

ptsInTriPolyh <- function(X,triPolyHedron)
{
     if(ncol(X)!=3) stop("X should be a 3-column matrix.")

   # Extract vertex and vertex difference vectors:
   
   n <- nrow(X)
   v1 <- triPolyHedron$v1
   v2 <- triPolyHedron$v2
   v3 <- triPolyHedron$v3

   Delta21 <- v2 - v1
   Delta32 <- v3 - v2

   # Obtain centroid:

   centroid <- triPolyhCentroid(triPolyHedron)

   # Remove redundant faces:
     
   condit <- Delta21==Delta32
   conditN <- cbind(as.numeric(condit[,1]),as.numeric(condit[,3]),
                    as.numeric(condit[,3]))
   conditS <- apply(conditN,1,sum)
   redFaceInds <- (1:length(conditS))[conditS>0]
   if (length(redFaceInds)>0)
   {
      v1 <- v1[-redFaceInds,] ; v2 <- v2[-redFaceInds,]
      v3 <- v3[-redFaceInds,] ; Delta21 <- Delta21[-redFaceInds,]
      Delta32 <- Delta32[-redFaceInds,]
   }

   # Obtain indicators for points inside
   # polyhedron using face-by-face comparison
   # with the centroid.
   
   A <- (Delta21[,c(2,3,1)]*Delta32[,c(3,1,2)]
        - Delta21[,c(3,1,2)]*Delta32[,c(2,3,1)])
   numFaces <- nrow(v1)
   indicX <- sign(tcrossprod(X,A) - 
                  matrix(rep(apply(A*v1,1,sum),n),n,numFaces,byrow=TRUE))
   indicC <- sign(apply(A*(t(centroid-t(v1))),1,sum))
   indics <- t(t(indicX)*indicC)

   return(apply(indics,1,sum)==numFaces)
}

############ End of ptsInTriPolyh #########

