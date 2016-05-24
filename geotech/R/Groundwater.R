##  GROUNDWATER
##  Kyle Elmy and Jim Kaklamanos
##  1 January 2016



########################################################################################
##  1. HYDRAULIC CONDUCTIVITY FROM CONSTANT HEAD TEST

##  Output:   Hydraulic condcutivity from the constant head test
##
##  Inputs:   V = volume of water collected
##            t = time of flow
##            h = head difference between inflow and outflow
##            L = length of soil sample
##            As = cross-sectional area of the soil sample
##            Ds = diameter of soil sample
##
##  Notes:    o  Either English or metric units can be used, but they must be consistent.
##            o  Either the area or the diameter of the soil sample need to be specified

kConstant <- function(V, t, h, L, As = NA, Ds = NA){

  ##  Compute flow rate
  Q <- V / t

  ##  Compute hydraulic gradient
  i <- h / L

  ##  Compute cross-sectional area
  if(is.na(As) == TRUE){
    As <- pi/4 * Ds^2
  }

  ##  Compute hydraulic conductivity
  k <- Q / (i*As)

  return(k)
}

##  Example code for this function (k in units of cm/s)
##  kConstant(V = 800, t = 100, h = 200, As = 40, L = 50)





########################################################################################
##  2. HYDRAULIC CONDUCTIVITY FROM FALLING HEAD TEST

##  Output:   Hydraulic condcutivity from the falling head test
##
##  Inputs:   h0 = head difference at beginning of test
##            hf = head difference at end of test (h0 > hf)
##            t = time of flow
##            L = length of soil sample
##            As = cross-sectional area of the soil sample
##            Ap = cross-sectional area of the standpipe
##            Ds = diameter of the soil sample
##            Dp = diameter of the standpipe
##
##  Notes:    o  Either English or metric units can be used, but they must be consistent.
##            o  Either the areas (As and Ap) OR the diameters (Ds and Dp) need to be specified.

kFalling <- function(h0, hf, t, L, As = NA, Ap = NA, Ds = NA, Dp = NA){

  ##  Cross-sectional areas
  if(is.na(As) == TRUE){
    As <- pi/4 * Ds^2
  }
  if(is.na(Ap) == TRUE){
    Ap <- pi/4 * Dp^2
  }

  ##  Hydraulic conductivity
  k <- ((Ap * L) / (As * t)) * log(h0 / hf)

  return(k)
}

##  Example code for this function (k in units of cm/s)
##  kFalling(h0 = 12, hf = 2, L = 10, Ds = 20, Dp = 2, t = 100)




########################################################################################
##  3. EQUIVALENT HORIZONTAL HYDRAULIC CONDUCTIVITY

##  Output:   Equivalent horizontal hydraulic conductivity (kx) for layered soil deposits
##
##  Inputs:   thk = vector of layer thicknesses
##            k = vector of layer hydraulic conductivities


kx <- function(thk, k){

  return(sum(k * thk) / sum(thk))
  
}




########################################################################################
##  4. EQUIVALENT VERTICAL HYDRAULIC CONDUCTIVITY

##  Output:   Equivalent vertical hydraulic conductivity (kz) for layered soil deposits
##
##  Inputs:   thk = vector of layer thicknesses
##            k = vector of layer hydraulic conductivities


kz <- function(thk, k){

  return(sum(thk) / sum(thk / k))
  
}




########################################################################################
##  5. FLOW RATE TO WELLS

##  Output:   Q = flow rate to well
##
##  Inputs:   k = hydraulic conductivity of aquifer
##            H = thickness of aquifer
##            h0 = initial total head in aquifer (before pumping)
##            hf = final total head in well casing (after pumping)
##            r0 = radius of influence
##            rw = radius of well
##
##  Notes:    o  Datum for h0 and hf is the bottom of the aquifer.
##            o  For unconfined aquifers, H > h0 (or specify as NA)
##            o  For confined aquifers, H >= hf.
##            o  For mixed aquifers (which start as confined prior to pumping and finish as
##               unconfined after pumping is complete), H < hf.

wellFlow <- function(k, H = NA, h0, hf, r0, rw){

  ##  Unconfined aquifers
  if(is.na(H) == TRUE || H > h0){
    Q <- pi * k * (h0^2 - hf^2) / log(r0 / rw)
  } else{
    if(is.na(H) == FALSE){

      ##  Confined aquifers
      if(H <= hf){
        Q <- 2 * pi * k * H * (h0 - hf) / log(r0 / rw)
      } else{

        ##  Mixed aquifers
        if(H > hf){
          Q <- pi * k * (2*H*h0 - H^2 - hf^2) / log(r0 / rw)
        }
      }
    }
  }
  return(Q)
}

##  Example:
##  wellFlow(k = 0.065, H = 10, h0 = 21, hf = 15, r0 = 20, rw = 2)




########################################################################################
##  6. WELL DRAWDOWN

##  Output:   A two-element list containing:
##            h = height of groundwater surface a distance r from the well
##            dd = h0 - h = drawdown of groundwater surface a distance r from the well
##
##  Inputs:   Q = flow rate into well
##            k = hydraulic conductivity of aquifer
##            H = saturated thickness of aquifer
##            h0 = initial total head in aquifer (before pumping)
##            r0 = radius of influence
##            rw = radius of well
##            r = radius of interest
##
##  Notes:    o  Datum for total heads are the bottom of the aquifer.
##            o  For unconfined aquifers, H = NA.
##            o  For confined aquifers, H >= hf.
##            o  For mixed aquifers (which start as confined prior to pumping and finish as
##               unconfined after pumping is complete), H < hf.

wellDrawdown <- function(Q, k, H = NA, h0, r0, rw, r){

  ##  Return h = h0 and dd = 0 if r > r0
  if(r >= r0){
    h <- h0
    dd <- 0
  } else{
  
    ##  Unconfined aquifers
    if(is.na(H) == TRUE || H > h0){
      h <- sqrt(h0^2 - (Q/(pi*k)) * log(r0 / r))
      dd <- h0 - h
    } else{
      if(is.na(H) == FALSE){
        
        ##  Determine drawdown at well
        hf <- sqrt(2*H*h0 - H^2 - (Q/(pi*k)) * log(r0 / rw))
        
        ##  Confined aquifers
        if(H <= hf){
          h <- h0 - (Q / (2*pi*k*H)) * log(r0 / r)
          dd <- h0 - h
        } else{
          
          ##  Mixed aquifers
          if(H > hf){
            h <- sqrt(2*H*h0 - H^2 - (Q/(pi*k)) * log(r0 / r))
            dd <- h0 - h
          }
        }
      }
    }
  }
  return(list(h = h, dd = dd))
}

##  Example:
##  wellDrawdown(Q = 14.5, k = 0.065, H = 10, h0 = 21, r0 = 20, r = 2, rw = 2)




########################################################################################
##  7. HYDRAULIC CONDUCTIVITY FROM PUMPING TESTS

##  Output:   k = hydraulic conductivity of aquifer
##
##  Inputs:   Q = flow rate into well
##            H = thickness of aquifer
##            h1 = total head in farthest observation well
##            h2 = total head in nearest observation well
##            r1 = radius from pumped well to farthest observation well
##            r2 = radius from pumped well to nearest observation well
##
##  Notes:    o  Datum for h0 and hf is the bottom of the aquifer.
##            o  For unconfined aquifers, H > h1 (or specify as NA)
##            o  For confined aquifers, H >= h2.
##            o  For mixed aquifers (which start as confined prior to pumping and finish as
##               unconfined after pumping is complete), H < h2.

kPump <- function(Q, H = NA, h1, h2, r1, r2){

  ##  Unconfined aquifers
  if(is.na(H) == TRUE || H > h1){
    k <- Q * log(r1 / r2) / (pi * (h1^2 - h2^2))
  } else{
    if(is.na(H) == FALSE){

      ##  Confined aquifers
      if(H <= h2){
        k <- Q * log(r1 / r2) / (2 * pi * H * (h1 - h2))
      } else{
        
        ##  Mixed aquifers
        if(H > h2){
          k <- Q * log(r1 / r2) / (pi*(2*H*h1 - H^2 - h2^2))
        }
      }
    }
  }
  return(k)
}

##  Example:
##  kPump(Q = 14.5, H = 10, h1 = 20, h2 = 16, r1 = 16, r2 = 8)



