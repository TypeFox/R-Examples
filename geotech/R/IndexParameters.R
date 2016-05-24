##  PHASE DIAGRAMS AND INDEX PARAMETERS
##  Kyle Elmy and Jim Kaklamanos
##  1 December 2015


########################################################################################
##  1. PHASE DIAGRAMS PLOT

##  Description:  This function plots phase diagrams from weights (or masses)
##                and volumes of a soil sample
##
##  Output:   Plot of  phase diagram
##
##  Inputs:   Ws = weight of solids
##            Ww = weight of water
##            Vs = volume of solids
##            Vw = volume of water
##            Va = volume of air
##            W.unit = measurement unit of weights
##            V.unit = measurement unit of volumes
##            mass = logical variable: TRUE for masses or FALSE for weights (default)
##
##  Notes:  o If any parameters are zero, please enter "0"; do not leave them blank.
##          o Plot is currently not to scale; may be edited in the future
##            to allow this functionality.

phase.plot <- function(Ws, Ww, Vs, Vw, Va, W.unit, V.unit, mass = FALSE){

  ##  Reference values
  x <- 0
  y <- 0

  ##  Total weights
  Wt <- Ws + Ww
  Vt <- Vs + Vw + Va
  Vv <- Vw + Va
  
  ##  Start plot
  par(mgp = c(2.8, 1, 0), las = 1)
  plot(x, y, type="l", yaxs = "i", xaxs = "i", xaxt="n", yaxt="n", col = "black",
       lwd = FALSE, xlab = "", ylab = "", xlim = c(0,80), ylim = c(0,90), bty = "n",
       main = "Phase Diagram", sub = "Note: Sketch is not to scale.")
  
  ##  Arrows for weight side
  arrows(x0 = 15, y0 = 10, x1 = 15, y1 = 30, length = 0.1, angle = 30, code = 3,
         col = par("fg"), lty = NULL, xpd = FALSE)
  arrows(x0 = 15, y0 = 30, x1 = 15, y1 = 50, length = 0.1, angle = 30, code = 3,
         col = par("fg"), lty = NULL, xpd = FALSE)
  arrows(x0 = 7, y0 = 10, x1 = 7, y1 = 70, length = 0.1, angle = 30, code = 3,
         col = par("fg"), lty = NULL, xpd = FALSE)

  ##  Arrows for volume side
  arrows(x0 = 45, y0 = 50, x1 = 45, y1 = 70, length = 0.1, angle = 30, code = 3,
         col = par("fg"), lty = NULL, xpd = FALSE)
  arrows(x0 = 45, y0 = 30, x1 = 45, y1 = 50, length = 0.1, angle = 30, code = 3,
         col = par("fg"), lty = NULL, xpd = FALSE)
  arrows(x0 = 45, y0 = 10, x1 = 45, y1 = 30, length = 0.1, angle = 30, code = 3,
         col = par("fg"), lty = NULL, xpd = FALSE)
  arrows(x0 = 60, y0 = 30, x1 = 60, y1 = 70, length = 0.1, angle = 30, code = 3,
         col = par("fg"), lty = NULL, xpd = FALSE)
  arrows(x0 = 70, y0 = 10, x1 = 70, y1 = 70, length = 0.1, angle = 30, code = 3,
         col = par("fg"), lty = NULL, xpd = FALSE)

  ##  Text along sides
  text(x = 12, y = 19, labels = bquote(W[s] == .(round(Ws, 3))), srt = 90, col = "black")
  text(x = 12, y = 39, labels = bquote(W[w] == .(round(Ww, 3))), srt = 90, col = "black")
  text(x = 5, y = 40, labels = bquote(W[t] == .(round(Wt, 3))), srt = 90, col = "black")
  text(x = 48, y = 19, labels = bquote(V[s] == .(round(Vs, 3))), srt = 90, col = "black")
  text(x = 48, y = 39, labels = bquote(V[w] == .(round(Vw, 3))), srt = 90, col = "black")
  text(x = 48, y = 59, labels = bquote(V[a] == .(round(Va, 3))), srt = 90, col = "black")
  text(x = 62, y = 45, labels = bquote(V[v] == .(round(Vv, 3))), srt = 90, col = "black")
  text(x = 72, y = 40, labels = bquote(V[t] == .(round(Vt, 3))), srt = 90, col = "black")
  text(x = 30, y = 20, labels = "SOLIDS", col= "brown")
  text(x = 30, y = 40, labels = "WATER", col = "blue")
  text(x = 30, y = 60, labels = "AIR", col = "grey20")
  
  ##  Text above each side
  if(mass == TRUE){
    text(x = 10, y = 80, labels = paste("Masses (", W.unit, ")", sep = ""))
  } else{
    if(mass == FALSE){
      text(x = 10, y = 80, labels = paste("Weights (", W.unit, ")", sep = ""))
    }
  }
  text(x = 60, y = 80, labels = paste("Volumes (", V.unit, ")", sep = ""))

  ##  Line for Limits
  abline(h = segments(x0 = 20, x1 = 20, y0 = 10, y1 = 70))
  abline(h = segments(x0 = 40, x1 = 40, y0 = 10, y1 = 70))
  abline(h = segments(x0 = 20, x1 = 40, y0 = 70, y1 = 70))
  abline(h = segments(x0 = 20, x1 = 40, y0 = 10, y1 = 10))
  abline(h = segments(x0 = 20, x1 = 40, y0 = 30, y1 = 30))
  abline(h = segments(x0 = 20, x1 = 40, y0 = 50, y1 = 50))
  abline(h = segments(x0 = 5, x1 = 15, y0 = 78, y1 = 78))
  abline(h = segments(x0 = 55, x1 = 65, y0 = 78, y1 = 78))
  
}

##  Example:
##  phase.plot(Ws = 75.8, Ww = 15.6, Vs = 0.45, Vw = 0.25, Va = 0.1, W.unit = "lb",
##             V.unit = "ft^3", mass = FALSE)




########################################################################################
##  2. WEIGHT-VOLUME PARAMETERS

##  Description:  This function calculates a comprehensive list of index parameters from
##                provided weights and volumes of a soil sample.
##
##  Output:   The output is a ten-element list containing:
##            w = water content (as decimal)
##            S = degree of saturation (as decimal)
##            e = void ratio (as decimal)
##            n = porosity (as decimal)
##            Gs = specific gravity
##            gammaT = total unit weight
##            gammaD = dry unit weight
##            gammaS = unit weight of solids
##            gammaW = unit weight of water
##            gammaB = buoyant unit weight
##
##  Inputs:   Ws = weight of solids
##            Ww = weight of water
##            Vs = volume of solids
##            Vw = volume of water
##            Va = volume of air

phase.params <- function(Ws, Ww, Vs, Vw, Va){

  ##  Total weights and volumes
  Wt <- Ws + Ww
  Vt <- Vs + Vw + Va
  Vv <- Vw + Va

  ##  Water content
  w <- Ww / Ws

  ##  Degree of saturation
  S <- Vw / Vv

  ##  Void ratio
  e <- Vv / Vs

  ##  Porosity
  n <- Vv / Vt

  ##  Total unit weight
  gammaT <- Wt / Vt

  ##  Dry unit weight
  gammaD <- Ws / Vt

  ##  Unit weight of solids
  gammaS <- Ws / Vs

  ##  Unit weight of water
  gammaW <- Ww / Vw
  
  ##  Specific gravity
  Gs <- Ws / (Vs * gammaW)

  ##  Buoyant unit weight
  gammaB <- gammaT - gammaW

  ##  Return
  return(list(w = w, S = S, e = e, n = n, Gs = Gs, gammaT = gammaT,
              gammaD = gammaD, gammaS = gammaS, gammaW = gammaW, gammaB = gammaB))
}

##  Example:
##  phase.params(Ws = 75.8, Ww = 15.6, Vs = 0.45, Vw = 0.25, Va = 0.1)




########################################################################################
##  3. WATER CONTENT

##  Description:  Calculates water content from lab results (i.e., measured weights
##                or masses)
##
##  Output:   Water content as a decimal
##
##  Inputs:    M1 = Mass (or weight) of can and wet soil, before oven dry
##             M2 = Mass (or weight) of can and dry soil, after drying in oven
##             Mc = Mass (or weight) of can
##
##  Notes:  Either mass or weight can be used, since the units cancel.

waterContent <- function(M1, M2, Mc){

  ##  Mass of water
  Mw <- M1 - M2
  
  ##  Mass of dry soil
  Ms <- M2 - Mc

  ##  Water content
  w <- Mw / Ms

  return(w)
}

##  Example:
##
##  waterContent(M1 = 20.68, M2 = 18.14, Mc = 8.20)




########################################################################################
##  4. RELATIVE DENSITY

##  Output:   Relative density as a decimal
##
##  Inputs:   e = void ratio
##            emax = maximum void ratio
##            emin = minimum void ratio

relDensity <- function(e, emax, emin){

  ## Numerator
  N <- emax - e

  ## Denominator
  D <- emax - emin

  ## Relative density
  Dr <- N / D

  return(Dr)
}

##  Example:
##
##  relDensity(e = 0.3, emax = 0.92, emin = 0.35)

