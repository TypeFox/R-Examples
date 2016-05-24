"keyfct.hn" <-
function(distance, key.scale)
#
# keyfct.hn
#
# Half-normal key function: exp (- (x^2/(2*sigma^2))
# The 2 is included in the denominator to preserve the 
# relationship with the standard deviation in a normal distribution.
#
# Arguments:
#
# distance   - perpendicular distance vector
# key.scale  - scale parameters
#
# Value:
#
# The routine returns a vector of probabilities that the observation
# were detected given they were at the specified distance and assuming that
# g`(0)=1 (ie a standard line transect detection function).
#
{
  exp( - (( distance/ (sqrt(2) * key.scale) )^2) )
}

