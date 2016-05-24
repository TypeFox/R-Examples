######## R function: triPolyhCentroid #########

# For determining the centroid of a triangular-faced
# polyhedron.

# Last changed: 02 DEC 2008

triPolyhCentroid <- function(triPolyHedron)
{
   v1 <- triPolyHedron$v1
   v2 <- triPolyHedron$v2
   v3 <- triPolyHedron$v3
   Rvecs <- (v1+v2+v3)/3

   sumsq <- function(x)
      return(sum(x^2))

   Areas <- sqrt(apply(extprod3d((v2-v1),(v3-v1)),1,sumsq))

   return(apply(Areas*Rvecs,2,sum)/sum(Areas))
}

############ End of triPolyhCentroid #########

