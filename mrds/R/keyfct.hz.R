# Hazard rate key function
# Function has the form 1 - exp (- (x/scale) ^ (-shape)))
#
# distance perpendicular distance vector
# key.scale vector of scale values
# key.shape vector of shape values
#
# return vector of probabilities that the observations were detected given they were at the specified distance and assuming that g(0)=1 (ie a standard line transect detection function).
# documented in ?distpdf
keyfct.hz <- function(distance, key.scale, key.shape){
  return(1 - exp( - (distance/key.scale)^( - key.shape)))
}
