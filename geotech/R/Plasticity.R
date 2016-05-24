##  PLASTICITY
##  Kyle Elmy and Jim Kaklamanos
##  1 December 2015



########################################################################################
##  1. LIQUID LIMIT
##
##  Description:  Calculate and plot a soil's liquid limit (LL) using the flow curve
##
##  Output:   Liquid limit (percent) and plot of soil's flow curve (optional)
##
##  Inputs:   N = vector of number of blows obtained from the liquid limit test
##            w = vector of water contents (PERCENT, not decimal)
##            draw = logical variable: TRUE if plot of flow curve is desired, FALSE to
##                   suppress the creation of the plot

LL <- function(N, w, draw){
  
  ##  Determine best-fit line for w, N data
  lmLL <- lm(w ~ log(N))

  ##  Calculate LL at N = 25
  a <- as.numeric(lmLL$coefficients[1])
  b <- as.numeric(lmLL$coefficients[2])
  LL <- a + b * log(25)
  
  ##  Create plot
  if(draw == TRUE){
    par(las = 1, cex = 1.2)
    plot(x = N, y = w, log = "x", xaxt = "n", col = "white",
         xlab = "Number of blows, N", ylab = "Water content, w (%)",
         main = "Flow Curve from Liquid Limit Test")
    logAxis(x = N)
    axis(side = 1, tick = FALSE, at = c(1, 20, 50, 100))
    points(x = N, y = w, cex = 1.3, lwd = 2)
    x.line = seq(from = min(N), to = max(N), length.out = 2)
    y.line = a + b*log(x.line)
    lines(x.line, y.line, col = "blue", lwd = 2)
    segments(x0 = 25, y0 = min(w)/2, x1 = 25, y1 = LL, lwd = 1, lty = 2, col = "blue")
    segments(x0 = min(N)/2, y0 = LL, x1 = 25, y1 = LL, lwd = 1, lty = 2, col = "blue")
    text(x = min(N), y = LL*1.05, adj = 0, labels = bquote("LL = "*.(round(LL, 0))), cex = 1.2)
  }

  ##  Return LL
  return(LL)
 
}

##  Example code for liquid limit analysis
##  LL(N = c(72, 37, 14), w = c(7, 15, 21), draw = TRUE)





########################################################################################
##  2. PLASTICITY INDEX
##
##  Description:  Calculates plasticity index from liquid limit and plastic limit
##
##  Output:  A two-element list composed of:
##           PI = Plasticity index (percent)
##           descr = qualitative description of soil based on Sowers (1979)
##
##  Inputs:  LL = Liquid Limit (percent)
##           PL = Plastic Limit (percent)
##
PI <- function(LL, PL){

  ##  Plasticity index
  PI <- LL-PL

  ##  Characteristics
  if(PI <= 3){
    descr <- "Nonplastic"
  } else if(PI > 3 && PI <= 15){
    descr <- "Slightly plastic"
  } else if(PI > 15 && PI <= 30){
    descr <- "Medium plastic"
  } else if(PI > 30){
    descr <- "Highly plastic"
  }
  
  return(list(PI = PI, descr = descr))
}

##  Example:
##
##  PI(LL=80, PL=30)




########################################################################################
##  3. LIQUIDITY INDEX
##
##  Description:  Calculates liquidity index from liquid limit, plastic limit,
##                and in-situ moisture content
##
##  Output:  A two-element list composed of:
##           LI = Liquidity index (percent)
##           descr = qualitative description of soil in its natural moisture state
##
##  Inputs:  w = in-situ moisture content (percent)
##           PL = Plastic Limit (percent)
##           LL = Liquid Limit (percent)

LI <- function(w, LL, PL){

  ## Liquidity index
  LI <- (w - PL) / (LL - PL)

  ##  Characteristics
  if(LI < 0){
    descr <- "Solid/Semisolid State"
  } else if(LI == 0){
    descr <- "At Plastic Limit"
  } else if(LI > 0 && LI < 1){
    descr <- "Plastic State"
  } else if(LI == 1){
    descr <- "At Liquid Limit"
  } else if(LI > 1){
    descr <- "Liquid State"
  }
            
  return(list(LI = LI, descr = descr))
}

##  Example:
##
##  LI(w = 55, PL = 20, LL = 50)




########################################################################################
##  4. CASAGRANDE'S PLASTICITY CHART
##
##  Output:  Plot of a soil's plasticity parameters on Casagrande's plasticity chart 
##
##  Inputs:  LL = Liquid Limit (percent)
##           PL = Plastic Limit (percent)
##           PI = Plasticity index (percent)
##
##  Notes:  Either PL or PI must be specified

plasticity.plot <- function(LL, PL = NA, PI = NA){

  ##  Calculate PI if not provided
  if(is.na(PI) == TRUE){
    PI <- LL - PL
  }

  ##  Expression for axes
  xLab <- "Liquid Limit, LL"
  yLab <- "Plasticity Index, PI"

  ##  Create plot
  par(mgp = c(2.8, 1, 0), las = 1)
  plot(LL, PI, type = "o", yaxs = "i", xaxs = "i", col = "blue", lwd = 2,
       xlim = c(0,100), ylim = c(0,70),
       xlab = xLab, ylab = yLab, main = "Casagrande's Plasticity Chart")
  
  ##  Label zones
  text(x = 70, y = 40, labels = "CH")
  text(x = 30, y = 12, labels = "CL")
  text(x = 70, y = 20, labels = "MH")
  text(x = 40, y = 5, labels = "ML")
  text(x = 22, y = 6, labels = "CL-ML")

  ##  Add lines
  ##  LL = 50 line
  segments(x0 = 50, x1 = 50, y0 = 0, y1 = 37.8)
  ##  U-line
  x <- seq(from = 16, to = 100, by = 0.5)
  U <- 0.9*(x - 8)
  lines(x, U, lty = 2)
  ##  A-line
  x <- seq(from = 25.7, to = 100.7, by = 1)
  A <- 0.7*(x - 20)
  lines(x, A, lty = 1)
  ##  Lines near initial portion of plot
  segments(x0 = 16, x1 = 30, y0 = 7, y1 = 7)
  segments(x0 = 16, x1 = 25.7, y0 = 4, y1 = 4)
  segments(x0 = 16, x1 = 16, y0 = 0, y1 = 7)
  ##  Labels
  text(x = 50, y = 42, srt = 45, labels = "U-Line")
  text(x = 58, y = 30, srt = 35, labels = "A-Line")

  ##  Plot point of interest
  points(LL, PI, col = "blue", pch = 16, cex = 1.2)

}

##  Example:
##  plasticity.plot(LL = 40, PL = 20)




########################################################################################
##  5. A-LINE
##
##  Output:  The plasticity index corresponding to the A-line on Casagrande's plasticity chart
##
##  Inputs:  LL = Liquid Limit (percent)
##
##  Notes:  Either PL or PI must be specified

A.line <- function(LL){
  if(LL <= 4){
    return(LL)
  } else if(LL > 4 && LL <= 25.5){
    return(4)
  } else if(LL > 25.5){
    return(0.73*(LL - 20))
  }
}
