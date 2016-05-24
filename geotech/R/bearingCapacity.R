##  BEARING CAPACITY
##  Kyle Elmy and Jim Kaklamanos
##  1 October 2015


########################################################################################
##  1. GROSS BEARING PRESSURE

##  Output:   Gross bearing pressure of footing (psf or kPa)
##
##  Inputs:   P = Vertical gross column load (lb or kN)
##            B = foundation width (ft or m)
##            L = foundation length (ft or m)
##            D = Depth of foundation (ft or m)
##            Dw = Depth of groundwater table below foundation base (ft or m)
##            metric = logical variable: TRUE (for metric units) or FALSE (for English units)
##            gammaW = unit weight of water (default = 62.4 pcf for English units; 9.81 kN/m^3 for metric units)
##            gammaC = unit weight of concrete (default = 150 pcf for English units; 23.6 kN/m^3 for metric units)
## 
##  Notes:   o  Either SI or English units can be used, but must stay consistent.
##           o  When specifying the length and width, L should be the longer of the two lengths.
##           o  For a continuous (strip) foundation, specify L = 1 and specify P as the load per unit length
##           o  When the groundwater table is deep or unknown, set Dw >= D.

bearingPressure <- function(P, B, L, D, Dw, metric, gammaW = NA, gammaC = NA){

  ##  Define gammaW and gammaC for metric units
  if(metric == TRUE){
    if(is.na(gammaW) == TRUE){
      gammaW <- 9.81
    }
    if(is.na(gammaC) == TRUE){
      gammaC <- 23.6
    }
  } else{
    ##  Define gammaW and gammaC for metric units
    if(metric == FALSE){
      if(is.na(gammaW) == TRUE){
        gammaW <- 62.4
      }
      if(is.na(gammaC) == TRUE){
        gammaC <- 150
      }
    }
  }
  

  ##  Ensure that pore water pressure is no less than zero
  u <- max(0, gammaW * (D - Dw))

  ##  Calculate bearing pressure
  q <- (P / (B*L)) + (gammaC * D) - u
  
  return(q)
}

##  Example code for this function:
##
##  BearingPressure(P = 1000, B = 2, L = 5, D = 6, Dw = 2, metric = FALSE)




########################################################################################
##  2. BEARING CAPACITY FACTORS

##  Output:   Bearing capacity factor (Nq, Nc, Ngamma) using either the
##            Terzaghi or Vesic methods
##
##  Inputs:   phi = friction angle (degrees)
##            case = "general" or "local" to indicate general or
##                   local shear failure ("general" is default)
##            method = "Terzaghi" or "Vesic" ("Terzaghi" is default)
##
##  Notes:    o  For local shear, the friction angle is reduced to a value equal to
##               atan(2/3 * tan(phi)).
##            o  Terzaghi's Ngamma is approximated by Coduto (2001).


##  ------------------------------------------------------
##  FUNCTION FOR Nq

Nq <- function(phi, case = "general", method = "Terzaghi"){

  ##  Input checking
  if(case != "general" && case != "local")
    stop("Case may only be 'general' or 'local'.")
  if(method != "Terzaghi" && method != "Vesic")
    stop("Case may only be 'Terzaghi' or 'Vesic'.")
  
  ##  Adjustment for local shear
  if(case == "local") phi <- atan(2/3 * tan(phi*pi/180))*180/pi

  ##  Angles
  phi.deg <- phi
  phi.rad <- phi * pi/180

  ##  Evaluate bearing capacity
  if(method == "Terzaghi"){
    num <- exp(2*pi*(3/4 - phi.deg/360)*tan(phi.rad))
    den <- 2*(cos((45 + phi.deg/2)*pi/180))^2
    Nq.value <- num/den
  } else{
    if(method == "Vesic"){
      Nq.value <- tan((45 + phi.deg/2)*pi/180)^2 * exp(pi * tan(phi.rad))
    }
  }
  return(Nq.value)
}

 
##  -----------------------------------------------------
##  FUNCTION FOR Nc

Nc <- function(phi, case = "general", method = "Terzaghi"){

  ##  Input checking
  if(case != "general" && case != "local")
    stop("Case may only be 'general' or 'local'.")
  if(method != "Terzaghi" && method != "Vesic")
    stop("Case may only be 'Terzaghi' or 'Vesic'.")
  
  ##  Adjustment for local shear
  if(case == "local") phi <- atan(2/3 * tan(phi*pi/180))*180/pi

  ##  Angles
  phi.deg <- phi
  phi.rad <- phi * pi/180

  ##  Evaluate bearing capacity
  if(method == "Terzaghi"){
    Nq.value <- Nq(phi = phi.deg, case = "general", method = "Terzaghi")
    if(phi.deg == 0) return(5.7)
    if(phi.deg > 0) return((Nq.value-1)/tan(phi.rad))
  } else{
    if(method == "Vesic"){
      Nq.value <- Nq(phi = phi.deg, case = "general", method = "Vesic")
      if(phi.deg == 0) return(5.14)
      if(phi.deg > 0) return((Nq.value-1)/tan(phi.rad))
    }
  }
}


##  -----------------------------------------------------
##  FUNCTION FOR Ngamma

Ngamma <- function(phi, case = "general", method = "Terzaghi"){

  ##  Input checking
  if(case != "general" && case != "local")
    stop("Case may only be 'general' or 'local'.")
  if(method != "Terzaghi" && method != "Vesic")
    stop("Case may only be 'Terzaghi' or 'Vesic'.")
  
  ##  Adjustment for local shear
  if(case == "local") phi <- atan(2/3 * tan(phi*pi/180))*180/pi

  ##  Angles
  phi.deg <- phi
  phi.rad <- phi * pi/180

  ##  Evaluate bearing capacity
  if(method == "Terzaghi"){
    Nq.value <- Nq(phi = phi.deg, case = "general", method = "Terzaghi")
    num <- 2*(Nq.value+1)*tan(phi.rad)
    den <- 1 + 0.4*sin(4*phi.rad)
    Ngamma.value <- num/den
  } else{
    if(method == "Vesic"){
      Nq.value <- Nq(phi = phi.deg, case = "general", method = "Vesic")
      Ngamma.value <- 2*(Nq.value + 1)*tan(phi.rad)
    }
  }
  return(Ngamma.value)
}




########################################################################################
##  3. ULTIMATE BEARING CAPACITY (Terzaghi)


##  Output:   Bearing capacity (q_ult) from Terzaghi's simple theory (psf or kPa)
##
##  Inputs:   phi = effective friction angle (deg)
##            c = effective cohesion (psf or kPa)
##            B = foundation width (ft or m), or foundation diameter for circular footings
##            L = foundation length (ft or m)
##            D = Depth of foundation (ft or m)
##            Dw = Depth of groundwater table below foundation base (ft or m)
##            gamma = unit weight of soil (pcf of kN/m^3)
##            gammaW = unit weight of water (default = 62.4 pcf for English units; 9.81 kN/m^3 for metric units)
##            case = "general" or "local" to indicate general or
##                   local shear failure ("general" is default)
##            shape = "square", "rectangle", "circle", "strip" (or "continuous")
##            metric = logical variable: TRUE (for metric units) or FALSE (for English units)
## 
##  Notes:   o  Either SI or English units can be used, but must stay consistent.
##           o  When specifying the length and width of rectangular foundations,
##              L should be the longer of the two lengths.
##           o  When the groundwater table is deep or unknown, set Dw >= D + B

bearingCapacity <- function(phi, c, B, L, D, Dw, gamma, gammaW = NA, metric,
                            case = "general", shape = "square"){
  
  ##  Input checking
  if(case != "general" && case != "local")
    stop("Case may only be 'general' or 'local'.")
  if(shape != "rectangle" && shape != "square" && shape != "circle" &&
     shape != "strip" && shape != "continuous"){
    stop("Shape may only be 'rectangle', 'square', 'circle', or 'strip' (or 'continuous').")
  }

  ##  Define gammaW
  if(metric == TRUE){
    if(is.na(gammaW) == TRUE){
      gammaW <- 9.81
    }
  } else{
    if(metric == FALSE){
      if(is.na(gammaW) == TRUE){
        gammaW <- 62.4
      }
    }
  }
  
  ##  Correct cohesion for local shear
  if(case == "local"){
    c <- 2/3 * c
  }
  
  ##  Obtain bearing capacity factors
  Nc.value <- Nc(phi, case = case, method = "Terzaghi")
  Nq.value <- Nq(phi, case = case, method = "Terzaghi")
  Ngamma.value <- Ngamma(phi, case = case, method = "Terzaghi")

  ##  Effective stress at footing base
  sigmaD <- gamma * D -  max(0, gammaW * (D - Dw))
  
  ##  Correct unit weight for groundwater conditions
  ##  Case 1:  Groundwater table above base of footing
  if(Dw <= D){
    gamma <- gamma - gammaW
  } else{
    ##  Case 2:  Groundwater table within depth B of the footing base
    if(Dw > D && Dw < (D + B)){
      gamma <- gamma - gammaW*(1 - (Dw - D)/B)
    } else{
      ##  Case 3:  Groundwater table deeper than D + B below footing base
      if(Dw > (D + B)){
        gamma <- gamma
      }
    }
  }

  ##  Evaluate bearing capacity for different foundation shapes
  if(shape == "continuous" || shape == "strip"){
    qult <- (c * Nc.value) + (sigmaD * Nq.value) + (0.5 * gamma * B * Ngamma.value)
  } else{
    if(shape == "square"){
      qult <- (1.3 * c * Nc.value) + (sigmaD * Nq.value) + (0.4 * gamma * B * Ngamma.value)
    } else{
      if(shape == "circle"){
        qult <- (1.3 * c * Nc.value) + (sigmaD * Nq.value) + (0.3 * gamma * B * Ngamma.value)
      } else{
        if(shape == "rectangle"){
          qult <- (c * Nc.value * (1 + 0.2*(B/L))) + (sigmaD * Nq.value * (1 + 0.2*(B/L))) +
            (0.5 * gamma * B * Ngamma.value * (1 - 0.3*(B/L)))
        }
      }
    }
  }
  return(qult)
}
