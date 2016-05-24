##  FUNCTIONS FOR SITE EXPLORATION
##  Kyle Elmy and Jim Kaklamanos
##  1 January 2016


########################################################################################
##  1. PLOT OF SOIL PROFILE
##
##  Output:  A plot of a soil profile
##
##  Input:   thk = vector of layer thicknesses (ft or m)
##           depth = vector of layer bottom depths (ft or m)
##           zw = depth of groundwater table (ft or m)
##           type = vector of soil types (character strings)
##           gamma = vector of unit weights (pcf or kN/m^3)
##           phi = vector of soil friction angles (deg)
##           C = vector of soil cohesion (psf or kPa)
##           title = desired title of plot (deafult: "Soil Profile")
##           metric = logical variable: TRUE (for metric units) or FALSE (for English units)

##  Notes:   o  Either layer thicknesses or depth to layer bottoms must be specified.
##           o  The only necessary variables are thk (or depth), zw, and metric.  All other
##              variables are optional.

soil.profile <- function(thk = NA, depth = NA, zw, type = NA, gamma = NA,
                         phi = NA, C = NA, title = "Soil Profile", metric){

  ##  Determine depths
  if(all(is.na(depth) == TRUE)){
    depth <- cumsum(thk)
  }

  ##  Axis labels
  if(metric == TRUE){
    ylab <- "Depth (m)"
    gamma.unit <- " kN/m^3"
    c.unit <- " kPa"    
  } else{
    ylab <- "Depth (ft)"
    gamma.unit <- " pcf"
    c.unit <- " psf"
  }

  ##  Plot
  par(mgp = c(2.8, 1, 0), las = 1)
  plot(x = c(0, 1), y = c(0, max(depth)), ylim = c(max(depth), 0),
       col = "white", type="o", yaxs = "i", xaxs = "i", xaxt = "n",
       xlab = "", ylab = ylab,
       lwd = FALSE, main = "Soil Profile")
  
  midpt <- c()

  for (i in 1:length(depth)){
    abline(h=depth[i])
  }
  for(i in 1:length(depth)){
    if(i == 1){
      midpt[1] <- depth[1]/2
    } else{
      midpt[i] <- depth[i-1]+((depth[i]-depth[i-1])/2)
    }
  }
  for(i in 1:length(depth)){
    if(i == 1) depth.min <- 0 else depth.min <- depth[i-1]
    text(x = 0.5, y = depth.min + 0.5*(midpt[i]-depth.min),
         labels = type[i])
    count <- .5
    if(is.na(gamma[i]) == FALSE){
      text(x = 0.5, y = depth.min + 0.5*(midpt[i]-depth.min),
           pos = c(1, 2), offset = c(count, 0),
           adj = 0, labels = bquote(gamma==.(gamma[i])*.(gamma.unit)))
      count <- count + 1
    }
    if(is.na(phi[i]) == FALSE){
      text(x = 0.5, y = depth.min + 0.5*(midpt[i]-depth.min),
           pos = 1, offset = count,
           adj = 0, labels = bquote(phi==.(phi[i])*degree))
      count <- count + 1
    }
    if(is.na(C[i]) == FALSE){
      text(x = 0.5, y = depth.min + 0.5*(midpt[i]-depth.min),
           pos = 1, offset = count,
           adj = 0, labels = bquote(c==.(C[i])*.(c.unit)))
      count <- count + 1
    }
  }
    
  ##  Line for GWT
  if(all(is.na(zw)) == FALSE){
    abline(h = zw, lwd = 1, lty = 1, col = "white")
    abline(h = zw, lwd = 1, lty = 3, col="blue")
    points(x = 0.1, y = zw/1.05, pch = 6, cex = 1.1,)
  }
}

##  EXAMPLE: 
##  soil.profile(depth = c(20, 40, 52, 60), zw = 20,
##               type = c("Dry Sand", "Saturated Sand", "Soft Clay", "Dense Gravel"),
##               gamma = c(110, 115, 120, 150), phi = c(30, 30, NA, 38),
##               C = c(NA, NA, 300, NA), metric = FALSE)
 




########################################################################################
##  2. SPT CORRECTION: N60
##
##  Output:  SPT blow count (N-value) corrected for field procedures
##
##  Input:   N = raw SPT N-value
##           Lr = rod length (ft or m)
##           E = hammer efficiency as a decimal (default: 0.60)
##           Db = borehole diameter (in or mm)
##           SS = logical variable: TRUE for standard sampler [default]; FALSE
##                for sampler without liner
##           metric = logical variable: TRUE (for metric units) or FALSE (for English units)

N60 <- function(N, Lr, Db, SS = TRUE, E = 0.60, metric){

  ##  Convert units
  if(metric == FALSE){
    Db <- Db * 25.4
    Lr <- Lr / 3.28084
  }
    
  ##  Borehole diameter factor
  if(Db <= 130){
    Cb <- 1
  } else{
    if(Db > 130 && Db <= 175){
      Cb <- 1.05
    } else{
      if(Db > 175){
        Cb <- 1.15
      }
    }
  }

  ##  Sampling method factor
  if(SS == TRUE){
    Cs <- 1
  } else{
    if(SS == FALSE){
      Cs <- 1.2
    }
  }

  ##  Rod length factor
  if(Lr < 4){
    Cr <- 0.75
  } else{
    if(Lr == 4){
      Cr <- 0.8
    } else{
      if(Lr > 4 && Lr < 6){
        Cr <- 0.85
      } else{
        if(Lr == 6){
          Cr <- 0.9
        } else{
          if(Lr > 6 && Lr < 10){
            Cr <- 0.95
          } else{
            if(Lr == 10){
              Cr <- 0.975
            } else{
              if(Lr > 10){
                Cr <- 1
              }
            }
          }
        }
      }
    }
  }
  
  ##  N60 Correction
  N60 <- N * Cb * Cs * Cr * (E / 0.60)

  ##  Return N60 Value
  return(N60)
}

##  Example code for this function:
##  N60(N = 11, Lr = 25, Db = 4, E = 0.50, SS = TRUE, metric = FALSE)





########################################################################################
##  3. SPT CORRECTION: N1,60
##
##  Output:  SPT blow count (N-value) corrected for field procedures and overburden pressure
##
##  Input:   N60 = SPT N-value corrected for field procedures
##           sigma = effective vertical stress at the depth of interest (psf or kPa)
##           metric = logical variable: TRUE (for metric units) or FALSE (for English units)

N160 <- function(N60, sigma, metric){

  ##  Determine atmospheric pressure
  if(metric == TRUE){
    pa <- 100
  } else{
    if(metric == FALSE){
      pa <- 2000
    }
  }
      
  ##  Determine N160 value
  N160 <- N60 * sqrt(pa / sigma)

  ## Return values
  return(N160)
}

##  Example code for this function:
##
##  In English units
##  N160(N60 = 8, sigma = 1500, metric = FALSE)
##
##  In SI units
##  N160(N60 = 8, sigma = 90, metric = TRUE) 


