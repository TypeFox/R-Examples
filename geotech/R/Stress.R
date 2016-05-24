##  STRESS
##  Kyle Elmy and Jim Kaklamanos
##  1 October 2015


########################################################################################
##  1. VERTICAL STRESS AT A POINT
##
##  Output:  A three-element list containing:
##           sigmaZ.eff = effective vertical stress
##           sigmaZ.total = total vertical stress
##           u = pore water pressure
##
##  Input:   gamma = vector of unit weights (pcf or kN/m^3)
##           thk = vector of layer thicknesses (ft or m)
##           depth = vector of layer bottom depths (ft or m)
##           zw = depth of groundwater table (ft or m)
##           zout = desired depth of output (ft or m): a single value
##           gammaW = unit weight of water (default = 62.4 pcf for English units;
##                                                   9.81 kN/m^3 for metric units)
##           metric = logical variable: TRUE (for metric units) or FALSE (for English units)
##
##  Note:    Either layer thicknesses or depths to layer bottoms must be specified.

sigmaZ <- function(gamma, thk = NA, depth = NA, zw, zout, gammaW = NA, metric){

  ##  Calculate layer thicknesses, if not provided
  if(all(is.na(thk) == TRUE)){
    if(length(depth) == 1){
      thk <- depth
    } else{
      if(length(depth) > 1){
        thk <- c(depth[1], diff(depth))
      }
    }
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

  ##  Case: zout = 0
  if(zout == 0){
    sigmaZ.total <- 0


  ##  Case: zout > 0
  } else{
  
    ##  Add ficticious zeroth layer
    thk <- c(0, thk)
    gamma <- c(0, gamma)
    
    ##  Depths of layer interfaces
    zint <- cumsum(thk)
    
    ##  Layer above point of interest
    Nlayers <- max(which(zint < zout))
    
    ##  Total stress
    sigmaZ.total <- sum(thk[1:Nlayers] * gamma[1:Nlayers]) + (zout - zint[Nlayers]) * gamma[Nlayers+1]
  }
  
  ##  Pore water pressure
  u <- gammaW * max(0, zout - zw)

  ##  Effective stress
  sigmaZ.eff <- sigmaZ.total - u
  
  return(list(sigmaZ.eff = sigmaZ.eff, sigmaZ.total = sigmaZ.total, u = u))
}

##  EXAMPLE:
##  sigmaZ(gamma = c(108, 116), depth = c(15, 40), zout = 18, zw = 15, metric = FALSE)




########################################################################################
##  2. VERTICAL STRESS PROFILE
##
##  Output:  A four-element list containing four vectors:
##           depth = depth
##           sigmaZ.eff = effective vertical stress
##           sigmaZ.total = total vertical stress
##           u = pore water pressure
##
##  Input:   gamma = vector of unit weights (pcf or kN/m^3)
##           thk = vector of layer thicknesses (ft or m)
##           depth = vector of layer bottom depths (ft or m)
##           zw = depth of groundwater table (ft or m)
##           zout = desired depths of output (ft or m): defaults to critical locations in
##                  the profile (the top and bottom of the profile, layer interfaces, and the
##                  groundwater table)  [recommended]
##           gammaW = unit weight of water (default = 62.4 pcf for English units;
##                                                   9.81 kN/m^3 for metric units)
##           metric = logical variable: TRUE (for metric units) or FALSE (for English units)

sigmaZ.profile <- function(gamma, thk = NA, depth = NA, zw, zout = NA, gammaW = NA, metric){
  
  ##  Determine critical depths
  if(all(is.na(depth) == TRUE)){
    depth <- cumsum(thk)
  }
  if(all(is.na(zout) == TRUE)){
    zout <- sort(unique(c(0, depth, zw)))
  }
  
  ##  Calculate stresses
  sigmaZ.eff <- vector(length = length(zout))
  sigmaZ.total <- vector(length = length(zout))
  u <- vector(length = length(zout))
  for(i in 1:length(zout)){
    temp <- sigmaZ(gamma = gamma, thk = thk, depth = depth, zw = zw,
                   zout = zout[i], gammaW = gammaW, metric = metric)
    sigmaZ.eff[i] <- temp$sigmaZ.eff
    sigmaZ.total[i] <- temp$sigmaZ.total
    u[i] <- temp$u
  }
 
  return(list(depth = zout, sigmaZ.eff = sigmaZ.eff, sigmaZ.total = sigmaZ.total, u = u))
}

##  EXAMPLE:
##  sigmaZ.profile(gamma = c(108, 116), depth = c(15, 40), zw = 15, metric = FALSE)




########################################################################################
##  3. VERTICAL STRESS PLOT
##
##  Output:  Plot of vertical stress versus depth
##
##  Input:   depth = vector of depths (ft or m)
##           sigmaZ.eff = vector of effective stresses (psf or kPa)
##           sigmaZ.total = vector of total stresses (psf or kPa)
##           u = vector of pore water pressures (psf or kPa)
##           metric = logical variable: TRUE (for metric units) or FALSE (for English units)
##
##  Notes:   o  If sigmaZ.total and u are left blank, the plot is only constructed for effective stress
##           o  Once constructed, additional profiles may be added to this plot
##              (for example, for induced stress or maximum past pressure)

sigmaZ.plot <- function(depth, sigmaZ.eff, sigmaZ.total = NA, u = NA, metric){
  
  ##  Expression for axes
  if(metric == TRUE){
    xLab <- "Stress (kPa)"
    yLab <- "Depth, z (m)"
  } else{
    if(metric == FALSE){
      xLab <- "Stress (psf)"
      yLab <- "Depth, z (ft)"
    }
  }

  ##  X axis limit
  if(all(is.na(sigmaZ.total)) == FALSE){
    xMax <- max(sigmaZ.total)
  } else{
    xMax <- max(sigmaZ.eff)
  }
  
  ##  Plot
  par(mgp = c(2.8, 1, 0), las = 1)
  plot(sigmaZ.eff, depth, type = "l", yaxs = "i", xaxt = "n", col = "black", lwd = 2,
       xlim = c(0, xMax), ylim = c(max(depth), 0),
       xlab = "", ylab = yLab)
  axis(3)
  mtext(xLab, side = 3, line = 2.5)
  if(all(is.na(sigmaZ.total)) == FALSE && all(is.na(u)) == FALSE){
    lines(u, depth, col = "blue", lwd = 1)
    lines(sigmaZ.total, depth, col = "gray70", lty = 1, lwd = 3)
    lines(sigmaZ.eff, depth, col = "black", lwd = 2)
  }
 
  ##  Line for gwt
  if(all(is.na(u)) == FALSE){
    Zgwt <- depth[max(which(u == 0))]
    abline(h = Zgwt, lwd = 1, lty = 3)
    points(x = 0.95*xMax, y = Zgwt*0.97, pch = 6, cex = 1.1)
  }
    
           
  ## Legend
  if(all(is.na(sigmaZ.total)) == FALSE && all(is.na(u)) == FALSE){
    legend(x = "topright", bty = "n", lty = c(1, 1, 1), lwd = c(3, 2, 1),
           col = c("gray70", "black", "blue"),
           legend = c("Total vertical stress", "Effective vertical stress",
             "Pore water pressure"))
  } else{
    legend(x = "topright", bty = "n", lty = c(1), lwd = c(3),
           col = c("black"), legend = c("Effective vertical stress"))
  }
}


##  EXAMPLE:
##  Site with constant unit weight = 100 pcf, GWT at 10 ft depth
##  temp <- sigmaZ.profile(gamma = rep(100, 3), depth = c(10, 20, 30), zw = 10, metric = FALSE)
##  depth <- temp$depth
##  sigmaTotal <- temp$sigmaZ.total
##  u <- temp$u
##  sigmaEff <- temp$sigmaZ.eff
##  sigmaZ.plot(depth = depth, sigmaZ.eff = sigmaEff,
##              metric = FALSE, sigmaZ.total = sigmaTotal, u = u)




########################################################################################
##  4. INDUCED STRESS DUE TO POINT LOAD
##
##  Output:  A six-element list containing the induced stresses due to a point load
##           following the Boussinesq method:
##           o  sigmaX = induced horizontal stress in the x direction
##           o  sigmaY = induced horizontal stress in the y direction
##           o  sigmaZ = induced vertical stress
##           o  tauZX = induced shear stress on the ZX plane
##           o  tauYX = induced shear stress on the YX plane
##           o  tauYZ = induced shear stress on the YZ plane
##
##  Input:   P = point load
##           x = horizontal distance from load in the X direction
##           y = horizontal distance from load in the Y direction
##           z = depth of interest
##           nu = Poisson's ratio

induced.point <- function(P, x, y, z, nu){

  ##  Compute Pythagorean distances
  R <- sqrt((x^2)+(y^2)+(z^2))
  r <- sqrt((x^2)+(y^2))

  ##  Factor
  f <- P/(2*pi)
  
  ##  -----------------------------------------------
  ##  Induced horizontal stress in X direction
  ##  First section
  S1X <- (3*(x^2)*z)/(R^5)
  ##  Second section
  S2X <- ((x^2)-(y^2))/((R*(r^2))*(R+z))
  ##  Third section
  S3X <- ((y^2)*z)/((R^3)*(r^2))
  ##  Combine second and third
  S4X <- (1-(2*nu))*(S2X+S3X)
  ##  Compute change in stress in the x direction
  sigmaX <- f*(S1X-S4X)

  ##  -----------------------------------------------
  ##  Induced horizontal stress in Y direction
  ##  First section
  S1Y <- ((3*(y^2)*z)/(R^5))
  ##  Second section
  S2Y <-((y^2)-(x^2))/((R*(r^2))*(R+z))
  ##  Third section
  S3Y <- ((x^2)*z)/((R^3)*(r^2))
  ##  Combine second and third
  S4Y <- (1-(2*nu))*(S2Y+S3Y)
  ##  Compute change in stress in the y direction
  sigmaY <- f*(S1Y-S4Y)

  ##  -----------------------------------------------
  ##  Induced vertical stress in Z direction
  sigmaZ <- (3*P*(z^3))/(2*pi*(R^5))

  ##  -----------------------------------------------
  ##  Induced shear stress in the ZX direction
  tauZX <- (3*P*(z^2)*x)/(2*pi*(R^5))

  ##  -----------------------------------------------
  ##  Induced shear stress in the YX direction
  ##  First section
  S1yx <- ((3*x*y*z)/(R^5))
  ##  Second Section
  S2yx <- (((2*R)+2)*x*y)/((R+(z^2))*(R^3))
  ##  Combine first and second section
  S3yx <- (1-(2*nu))*S2yx
  ##  Compute TauYX
  tauYX <- f*(S1yx-S3yx)
  
  ##  -----------------------------------------------
  ##  Induced shear stress in the YZ direction
  tauYZ <- (3*P*(z^2)*y)/(2*pi*(R^5))

  ##  Return elements
  return(list(sigmaX = sigmaX, sigmaY = sigmaY, sigmaZ = sigmaZ,
              tauZX = tauZX, tauYX = tauYX, tauYZ = tauYZ))
}

##  EXAMPLE:
##  induced.point(P = 100000, x = 5, y = 2, z = 6, nu = 0.35)




########################################################################################
##  5. INDUCED STRESS DUE TO POINT LOAD: PROFILE
##
##  Output:  A seven-element list containing seven vectors displaying the depth variation of
##           induced stresses due to a point load following the Boussinesq method:
##           o  depth = vector of depths
##           o  sigmaX = induced horizontal stress in the x direction
##           o  sigmaY = induced horizontal stress in the y direction
##           o  sigmaZ = induced vertical stress
##           o  tauZX = induced shear stress on the ZX plane
##           o  tauYX = induced shear stress on the YX plane
##           o  tauYZ = induced shear stress on the YZ plane
##
##  Input:   P = point load
##           x = horizontal distance from load in the X direction
##           y = horizontal distance from load in the Y direction
##           z = vector of depths of interest (default: 1-ft or 1-m increments,
##               to a maximum depth of 50 ft or 50 m)
##           nu = Poisson's ratio

induced.point.profile <- function(P, x, y, z = NA, nu){

  ##  Specify default values for depths
  if(all(is.na(z) == TRUE)){
    z <- seq(from = 0, to = 50, by = 1)
  }
   
  ##  Calculate stresses
  sigmaX <- vector(length = length(z))
  sigmaY <- vector(length = length(z))
  sigmaZ <- vector(length = length(z))
  tauZX <- vector(length = length(z))
  tauYX <- vector(length = length(z))
  tauYZ <- vector(length = length(z))  
  for(i in 1:length(z)){
    temp <- induced.point(P = P, x = x, y = y, z = z[i], nu = nu)
    sigmaX[i] <- temp$sigmaX
    sigmaY[i] <- temp$sigmaY
    sigmaZ[i] <- temp$sigmaZ
    tauZX[i] <- temp$tauZX
    tauYX[i] <- temp$tauYX
    tauYZ[i] <- temp$tauYZ
  }

  ##  Return elements
  return(list(depth = z, sigmaX = sigmaX, sigmaY = sigmaY, sigmaZ = sigmaZ,
              tauZX = tauZX, tauYX = tauYX, tauYZ = tauYZ))
}

##  EXAMPLE:
##  induced.point.profile(P = 100000, x = 5, y = 2, nu = 0.35)






########################################################################################
##  6. INDUCED STRESS DUE TO AREA LOAD
##
##  Output:  Induced vertical stress at the center of the loaded area
##
##  Input:   z = depth of interest
##           q = applied pressure at ground surface
##           B = width of loaded area
##           L = length of loaded area (rectangular foundations only)
##           shape = shape of loaded area (a string containing
##                   "circle", "square", "strip", or "rectangle")
##
##  Notes:   This function currently uses the approximate method of Poulos and Davis (1974).
##           More advanced formulations are expected in future versions of this package.

induced.area <- function(z, q, B, L = NA, shape){

  ##  Check input for shape
  if(shape != "circle" && shape != "square" && shape != "strip" && shape != "rectangle"){
    stop("Shape must be either 'circle', 'square', 'strip', or 'rectangle'.")
  }

  ## Compute the influence factor for circular loaded areas
  IC <- 1-(1/(1+(B/(2*z))^(2))^(1.5))

  ## Compute the influence factor for square loaded areas
  IS <- 1-(1/(1+(B/(2*z))^(2))^(1.76))
  
  ## Compute the influence factor for square loaded areas
  Ist <- 1-(1/(1+(B/(2*z))^(1.38))^(2.60))
            
  ## Compute the influence factor for rectangular loaded areas of width B and length L
  if(is.na(L) == FALSE){
    EXPI <- 1.38+((0.62*B)/L)
    EXPE <- 2.60-((0.84*B)/L)
    IR <- 1-(1/(1+((B/(2*z))^(EXPI)))^(EXPE))
  }

  ##  Induced stress
  if(shape == "circle"){
    I <- IC
  } else{
    if(shape == "square"){
      I <- IS
    } else{
      if(shape == "strip"){
        I <- Ist
      } else{
        if(shape == "rectangle"){
          I <- IR
        }
      }
    }
  }
  return(I * q)
}

##  EXAMPLE:
##  induced.area(z = 10, q = 1000, B = 3, shape = "square")




########################################################################################
##  7. INDUCED STRESS DUE TO AREA LOAD: PROFILE
##
##  Output:  A two-element list:
##           depth = vector of depths
##           sigmaZ = vector of induced vertical stresses below center of the loaded area
##
##  Input:   z = vector of depths (default: 1-ft or 1-m increments,
##               to a maximum depth of 50 ft or 50 m)
##           q = applied pressure at ground surface
##           B = width of loaded area
##           L = length of loaded area (rectangular foundations only)
##           shape = shape of loaded area (a string containing
##                   "circle", "square", "strip", or "rectangle")
##
##  Notes:   This function currently uses the approximate method of Poulos and Davis (1974).
##           More advanced formulations are expected in future versions of this package.

induced.area.profile <- function(z = NA, q, B, L = NA, shape){
 
  ##  Specify default values for depths
  if(all(is.na(z) == TRUE)){
    z <- seq(from = 0, to = 50, by = 1)
  }
  
  ##  Calculate stresses
  sigmaZ <- vector(length = length(z))
  for(i in 1:length(z)){
    sigmaZ[i] <- induced.area(z = z[i], q = q, B = B, L = L, shape = shape)
  }
  
  ##  Return values
  return(list(depth = z, sigmaZ = sigmaZ))
}

##  EXAMPLE:
##  induced.area.profile(q = 1000, B = 3, shape = "square")

