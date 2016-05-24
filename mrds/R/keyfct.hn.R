# Half-normal key function
# Function has the form exp (- (x^2/(2*sigma^2))
# The 2 is included in the denominator to preserve the
# relationship with the standard deviation in a normal distribution.
#
# distance perpendicular distance vector
# key.scale vector of scale values
#
# return a vector of probabilities that the observation were detected given they were at the specified distance and assuming that g`(0)=1 (ie a standard line transect detection function)
# documented in ?distpdf
keyfct.hn <- function(distance, key.scale){
  exp( - (( distance/ (sqrt(2) * key.scale) )^2) )
}
