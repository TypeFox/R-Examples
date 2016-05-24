##  SLOPE STABILITY
##  Jim Kaklamanos and Kyle Elmy
##  1 January 2016




########################################################################################
##  1. SLOPE STABILITY (INFINITE SLOPE ANALYSIS)

##  Output:  Factor of safety against shear failure on slopes using infinite slope analyses
##
##  Input:   Soil parameters:
##           o  c = soil cohesion
##           o  phi = soil friction angle (degrees)
##           o  gamma = soil unit weight
##           o  gammaW = unit weight of water (default = 62.4 pcf for English units;
##                                                       9.81 kN/m^3 for metric units)
##           Slope geometry:
##           o  alpha = slope angle (degrees)
##           o  D = depth to failure plane
##           o  zw = distance of groundwater table above failure plane
##                   (use 0 for a dry slope and D for a submerged slope with parallel seepage)
##
##           Units:
##           o  metric = logical variable: TRUE (for metric units: kN/m^3)
##                       or FALSE (for English units: pcf)
##                       [this is needed if gammaW is unspecified]
##           
##
##  Notes:   Assumptions of infinite slope analyses include (Coduto et al., 2011):
##           1. The slope face is planar and of infinite extent.
##           2. The failure surface is parallel to the slope face.
##           3. Vertical columns of equal dimensions through the slope are identical.
##

FSinf <- function(c, phi, gamma, gammaW = NA, alpha, D, zw, metric){

  
  ##  Define gammaW (if unspecified)
  if(is.na(gammaW) == TRUE){
    if(metric == TRUE){
      gammaW <- 9.81
    } else{
      if(metric == FALSE){
        gammaW <- 62.4
      }
    }
  }

  ##  Convert angles to radians
  alpha <- alpha * pi/180
  phi <- phi * pi/180

  ##  Numerator
  num <- c + (gamma*D - gammaW*zw) * (cos(alpha)^2) * tan(phi)

  ##  Demonimator
  den <- gamma * D * sin(alpha) * cos(alpha)

  ##  Factor of safety
  FS <- num / den
  
  return(FS)
}




########################################################################################
##  2. SLOPE STABILITY (PLANAR FAILURE ANALYSIS)

##  Output:  Factor of safety against shear failure on slopes with a planar failure surface
##
##  Input:   Soil parameters:
##           o  c = soil cohesion
##           o  phi = soil friction angle
##
##           Slope geometry:
##           o  alpha = angle of failure plane (degrees)
##           o  L = length of failure plane
##           o  W = weight of slope per unit width
##           o  u = average pressure head on the failure plane
##
##  Notes:   Either English or metric units can be used, but they must be consistent.

FSplanar <- function(c, phi, alpha, L, W, u){

  ##  Convert angles to radians
  alpha <- alpha * pi/180
  phi <- phi * pi/180

  ##  Numerator
  num <- c*L + (W*cos(alpha) - u*L) * tan(phi)

  ##  Demonimator
  den <- W * sin(alpha)

  ##  Factor of safety
  FS <- num / den
  
  return(FS)
}

