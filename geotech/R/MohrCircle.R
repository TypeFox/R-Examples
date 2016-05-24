##  MOHR CIRCLE ANALYSES
##  Kyle Elmy and Jim Kaklamanos
##  1 October 2015


########################################################################################
##  1. STRESS TRANSFORMATION

##  Output:   A two-element list containing:
##            sigma = normal stress on an inclined plane
##            tau = shear stress on an inclined plane
##
##  Inputs:   theta = angle of inclination (degrees)
##            sigmaX = normal stress acting in the horizontal direction
##            sigmaZ = normal stress acting in the vertical direction
##            tauXZ = shear stress acting on the same plane as sigmaX
##            sigma1 = major principal stress
##            sigma3 = minor principal stress

##  Notes:    o  In addition to theta, One of the following two sets of data must be entered:
##               1.  sigmaX, sigmaZ, tauXZ
##               2.  sigma1, sigma3
##            o  If theta is entered in conjunction with sigmaX, sigmaZ, and tauXZ, it is interpreted
##               as the angle of inclination above the horizontal.  If theta is entered in conjunction
##               with the principal stresses, then it is interpreted as the angle of inclination above the
##               major principal plane.

sigmaTrans <- function(theta, sigmaX = NA, sigmaZ = NA, tauXZ = NA, sigma1 = NA, sigma3 = NA){

  ##  Convert angle to radians
  theta <- theta * pi/180

  ##  Calculate normal and shear stresses from sigmaX, sigmaZ, tauXZ
  if(is.na(sigmaX) == FALSE && is.na(sigmaZ) == FALSE && is.na(tauXZ) == FALSE){
    sigma <- (sigmaZ + sigmaX)/2 + (sigmaZ - sigmaX)/2 * cos(2*theta) + tauXZ * sin(2*theta)
    tau <- (sigmaZ - sigmaX)/2 * sin(2*theta) - tauXZ * cos(2*theta)

    ##  Calculate normal and shear stresses from principal stresses
  } else{
    if(is.na(sigma1) == FALSE && is.na(sigma3) == FALSE){
      sigma <- (sigma1 + sigma3)/2 + (sigma1 - sigma3)/2 * cos(2*theta)
      tau <- (sigma1 - sigma3)/2 * sin(2*theta)
    }
  }

  return(list(sigma = sigma, tau = tau))
}

##  Example:
##  sigmaTrans(sigmaX = 80, sigmaZ = 120, tauXZ = 20, theta = 78)





########################################################################################
##  2. MOHR CIRCLE:  CALCULATE

##  Output:   A five-element list containing:
##            C = center of Mohr circle
##            R = radius of Mohr circle    
##            sigma = vector of normal stresses for Mohr circle
##            tau = vector of shear stresses for Mohr circle
##            theta = vector of angles (deg)
##
##  Inputs:   sigmaX = normal stress acting in the horizontal direction
##            sigmaZ = normal stress acting in the vertical direction
##            tauXZ = shear stress acting on the same plane as sigmaX
##            sigma1 = major principal stress
##            sigma3 = minor principal stress
##            theta = vector of angles (degrees); defaults to 0-180 in increments of 1

##  Notes:    o  One of the following two sets of data must be entered:
##               1.  sigmaX, sigmaZ, tauXZ
##               2.  sigma1, sigma3
##            o  If theta is entered in conjunction with sigmaX, sigmaZ, and tauXZ, it is interpreted
##               as the angle of inclination above the horizontal.  If theta is entered in conjunction
##               with the principal stresses, then it is interpreted as the angle of inclination above the
##               major principal plane.

MohrCircle.calc <- function(sigmaX = NA, sigmaZ = NA, tauXZ = NA, sigma1 = NA, sigma3 = NA,
                            theta = seq(from = 0, to = 180, by = 1)){
  
  ##  Calculate normal and shear stresses
  stress.vec <- sapply(X = theta, FUN = sigmaTrans, sigmaX = sigmaX, sigmaZ = sigmaZ,
                       tauXZ = tauXZ, sigma1 = sigma1, sigma3 = sigma3)
  sigma <- as.numeric(stress.vec[1,])
  tau <- as.numeric(stress.vec[2,])

  ##  Center and radius
  ##  Calculate from sigmaX, sigmaZ, tauXZ
  if(is.na(sigmaX) == FALSE && is.na(sigmaZ) == FALSE && is.na(tauXZ) == FALSE){
    C <- (sigmaX + sigmaZ) / 2
    R <- sqrt((((sigmaZ - sigmaX) /2) ^ 2) + ((tauXZ)^2))

    ##  Calculate from principal stresses
  } else{
    if(is.na(sigma1) == FALSE && is.na(sigma3) == FALSE){
      C <- (sigma1 + sigma3) / 2
      R <- (sigma1 - sigma3) / 2
    }
  }
  
  return(list(C = C, R = R, sigma = sigma, tau = tau, theta = theta))
}
  
##  Example
##  MohrCircle.calc(sigmaX = 80, sigmaZ = 120, tauXZ = 20)






########################################################################################
##  3. MOHR CIRCLE:  PLOT

##  Output:   Plot of Mohr Circle
##
##  Inputs:   sigmaX = normal stress acting in the horizontal direction
##            sigmaZ = normal stress acting in the vertical direction
##            tauXZ = shear stress acting on the same plane as sigmaX
##            sigma1 = major principal stress
##            sigma3 = minor principal stress
##            metric = logical variable: TRUE for metric units (kPa), and
##                     FALSE for English units
##
##  Notes:    o  One of the following two sets of data must be entered:
##               1.  sigmaX, sigmaZ, tauXZ
##               2.  sigma1, sigma3

MohrCircle.plot <- function(sigmaX = NA, sigmaZ = NA, tauXZ = NA, sigma1 = NA, sigma3 = NA,
                            metric = TRUE){
  
  ##  Calculate normal and shear stresses
  theta = seq(from = 0, to = 180, by = 1)
  stress.vec <- sapply(X = theta, FUN = sigmaTrans, sigmaX = sigmaX, sigmaZ = sigmaZ,
                       tauXZ = tauXZ, sigma1 = sigma1, sigma3 = sigma3)
  sigma <- as.numeric(stress.vec[1,])
  tau <- as.numeric(stress.vec[2,])
  
  ##  Expression for axes
  if(metric == TRUE){
    xLab <- expression("Normal Stress, "*sigma*" (kPa)")
    yLab <- expression("Shear Stress, "*tau*" (kPa)")
  } else{
    if(metric == FALSE){
      xLab <- expression("Normal Stress, "*sigma*" (psf)")
      yLab <- expression("Shear Stress, "*tau*" (psf)")
    }
  }
  
  par(mgp = c(2.8, 1, 0), las = 1)
  plot(sigma, tau, col = "black", type = "l", xlab = xLab, ylab = yLab, xaxs = "i",
       xlim = c(min(0, min(sigma)), max(sigma)*1.1), main = "Mohr Circle Plot", asp = 1)
  abline(h = 0)
}
  
##  Example
##  MohrCircle.plot(sigmaX = 80, sigmaZ = 120, tauXZ = 20, metric = FALSE)





########################################################################################
##  4. PRINCIPAL STRESSES

##  Output:   A four-element list containing:
##            sigma1 = magnitude of major principal stress
##            sigma3 = magnitude of minor principal stress
##            theta1 = direction of major principal stress (deg)
##            theta3 = direction of minor principal stress (deg)
##
##  Inputs:   sigmaX = normal stress acting in the horizontal direction
##            sigmaZ = normal stress acting in the vertical direction
##            tauXZ = shear stress acting on the same plane as sigmaX

sigma13 <- function(sigmaX, sigmaZ, tauXZ){

  ##  Major principal stress
  sigma1 <- ((sigmaX + sigmaZ) / 2) + (sqrt((((sigmaZ - sigmaX) /2) ^ 2) + ((tauXZ)^2)))
  ##  Angle of major principal plane
  x <- (1 / (1+((2*tauXZ)/(sigmaZ - sigmaX))^2))
  theta1 <- 1/2*(acos(sqrt(x)))*(180/pi)

  ##  Minor principal stress
  sigma3 <- ((sigmaX + sigmaZ) / 2) - (sqrt((((sigmaZ - sigmaX) /2) ^ 2) + ((tauXZ)^2)))
  ##  Angle of minor principal plane
  theta3 <- theta1 + 90

  ##  Return values
  return(list(sigma1 = sigma1, sigma3 = sigma3, theta1 = theta1, theta3 = theta3))
}

##  Example:
##  sigma13(sigmaX = 80, sigmaZ = 120, tauXZ = 20)





########################################################################################
##  5. MAXIMUM IN-PLANE SHEAR STRESS

##  Output:   A two-element list containing:
##            tauMax = maximum in-plane shear stress
##            theta = angle of maximum in-plane shear stress
##
##  Inputs:   sigmaX = normal stress acting in the horizontal direction
##            sigmaZ = normal stress acting in the vertical direction
##            tauXZ = shear stress acting on the same plane as sigmaX

tauMax <- function(sigmaX, sigmaZ, tauXZ){

  ##  Maximum in-plane shear stress
  tauMax <- sqrt((((sigmaZ - sigmaX) /2) ^ 2) + ((tauXZ)^2))
  ##  Angle of maximum in-plane shear stress
  x <- (1 / (1+((2*tauXZ)/(sigmaZ - sigmaX))^2))
  theta1 <- 1/2*(acos(sqrt(x)))*(180/pi)
  theta <- theta1 + 45

  ##  Return values
  return(list(tauMax = tauMax, theta = theta))
}

##  Example:
##  tauMax(sigmaX = 80, sigmaZ = 120, tauXZ = 20)
