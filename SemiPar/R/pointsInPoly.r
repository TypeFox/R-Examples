########## R-function: pointsInPoly ##########

# For determining those points are inside
# a polygon specified by "vertices"

# Last changed: 05 JUN 2009

pointsInPoly <- function(point.coords, vertices)
{
   if (vertices[1,1]!=vertices[nrow(vertices),1]) 
      vertices <- rbind(vertices,vertices[,1])
   if (vertices[1,2]!=vertices[nrow(vertices),2]) 
      vertices <- rbind(vertices,vertices[,1])

   n <- nrow(vertices) - 1

   # Translate to coords centered at each point

   x <-  - outer(point.coords[,1],vertices[1:n,1],"-")
   y <-  - outer(point.coords[,2],vertices[1:n,2],"-") 

   # Use tan(theta) = (v1 x v2) / (v1 * v2)

   i2 <- c(2:n, 1)

   arg1 <- x * y[, i2] - y * x[, i2]
   arg2 <-  x * x[, i2] + y * y[, i2]

   theta <- atan2(arg1,arg2)

   # Sum angles about each point

   theta <- abs(theta %*% rep(1, n))

   return(ifelse(abs(theta - 2*pi) < 1e-06, TRUE, FALSE))
}

########## End of pointsInPoly ##########
