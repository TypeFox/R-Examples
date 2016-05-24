################# R function: areaPolygon ################# 

# For obtaining the area inside a Polygon in dimension d=2 

# Last changed: 04 SEP 2008

areaPolygon <- function(polyVerts) {

   lenVal <- length(polyVerts$x)

   # Check equality of first and last vertices:
   if ( (polyVerts$x[1]!=polyVerts$x[lenVal]) | 
       (polyVerts$y[1]!=polyVerts$y[lenVal]) )
      stop("Start and end vertices must be equal.")

   # Obtain area.
   xVec <- polyVerts$x
   yVec <- polyVerts$y
		  
   return(abs(sum(xVec[1:(lenVal-1)]*yVec[2:lenVal]
          -xVec[2:lenVal]*yVec[1:(lenVal-1)])/2))

}

################# End of areaPolygon ################# 
