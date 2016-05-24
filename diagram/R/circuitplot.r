## =============================================================================
## Functions for drawing circuit elements
## =============================================================================
tpos <- function(mid,
                 pos,
                 dtext = 0.2) {
  pp  <- mid                
  dtext.x <- dtext.y <- dtext    # change for nonsymmetric plots
  adj <- c(0.5, 0.5)
  if (pos == 0) {
    pp <- mid 
    adj <- c(0.5, 0.5)
  } else if (pos == 1) {
    pp <- mid - c(0, dtext)
    adj <- c(0.5, 0)
  } else if (pos == 1.5) {
    pp <- mid - c(dtext, dtext/2)
    adj <- c(0, 1)
  } else if (pos == 2) {
    pp <- mid - c(dtext.x, 0)
    adj <- c(0, 0.5)
  }else if (pos == 2.5) {
    pp <- mid + c(-dtext.x, dtext.y/2)
    adj <- c(0, 0)
  } else if (pos == 3) {
    pp <- mid + c(0, dtext.y)
    adj <- c(0.5, 1)
  }else if (pos == 3.5) {
    pp <- mid + c(dtext.x, dtext.y/2)
    adj <- c(1, 0)
  } else if (pos == 4) {
    pp <- mid + c(dtext.x, 0)
    adj <- c(1, 0.5)
  }else if (pos == 4.5) {
    pp <- mid + c(dtext.x, -dtext.y/2)
    adj <- c(1, 1)
  }
  
  list (x = pp[1], y = pp[2], adj = adj)
} # end of tpos 

## -----------------------------------------------------------------------------

en.Resistor <- function (mid,
                         width = 0.05,  # the width of the resistor
                         length = 0.1,  # the length of the resistor
                         lab = NULL,    # the label
                         pos = 0,       # relative position (shift)
                         dtext = 0.,    # shift in x-y- coordinates for text
                         vert = TRUE,
                         ...) {
 if ( vert)
   rect(mid[1] - width/2, mid[2] - length/2,
        mid[1] + width/2, mid[2] + length/2, col = "white")
 else
   rect(mid[1] - length/2, mid[2] - width/2,
        mid[1] + length/2, mid[2] + width/2, col = "white")

  if (! is.null(lab))  {
    pp <- tpos(mid, pos, dtext = dtext)
    text(pp$x, pp$y, lab = lab, adj = pp$adj, ...)
  }
}
## -----------------------------------------------------------------------------

en.Capacitator <- function (mid,
                         width = 0.025,
                         length = 0.1,
                         lab = NULL,
                         pos = 2.5,
                         dtext = 0.04, 
                         vert = TRUE,
                         ...) {
  if (! vert) {
    rect(mid[1] - width/2, mid[2] - length/2,
         mid[1] + width/2, mid[2] + length/2,
         col = "white", border = NA)
    segments(mid[1] - width/2, mid[2] - length/2,
             mid[1] - width/2, mid[2] + length/2,
             lwd = 2)
    segments(mid[1] + width/2, mid[2] - length/2,
             mid[1] + width/2, mid[2] + length/2,
             lwd = 2)
  } else {
    rect(mid[1] - length/2, mid[2] - width/2,
         mid[1] + length/2, mid[2] + width/2,
         col = "white", border = NA)
    segments(mid[1] - length/2, mid[2] - width/2,
             mid[1] + length/2, mid[2] - width/2,
             lwd = 2)
    segments(mid[1] - length/2, mid[2] + width/2,
             mid[1] + length/2, mid[2] + width/2,
             lwd = 2)
  }
  if (! is.null(lab))  {
    pp <- tpos(mid, pos, dtext = dtext)
    text(pp$x, pp$y, lab = lab, adj = pp$adj, ...)
  }   
}

## -----------------------------------------------------------------------------

en.Node <- function (mid,
                  cex = 1,
                  lab = NULL,
                  pos = 2.5,
                  dtext = 0.025,
                  ...) {
  points(mid[1], mid[2], cex = cex, pch = 16, ...)
  if (! is.null(lab))  {
    pp <- tpos(mid, pos, dtext = dtext)
    text(pp$x, pp$y, lab = lab, adj = pp$adj, ...)
  }   
}

## -----------------------------------------------------------------------------

en.Amplifier <- function (mid,
                       r = 0.05,
                       lab = NULL,
                       pos = 0,
                       dtext = 0,
                       ...) {
  plotcircle(mid = mid, r = r)
  segments(mid[1], mid[2] - r/2,
           mid[1], mid[2] + r/2)
  if (! is.null(lab))  {
    pp <- tpos(mid, pos, dtext = dtext)
    text(pp$x, pp$y, lab, cex = pp$adj, ...)
  }   
}  

## -----------------------------------------------------------------------------

en.Transistor <- function (mid, gate, drain, source,
                       r = 0.05,
                       lab = NULL,
                       pos = 0,
                       dtext = 0,
                       ...) {
  segments(gate[1], gate[2], mid[1], mid[2])
  segments(drain[1], drain[2], mid[1], mid[2])
  segments(source[1], source[2], mid[1], mid[2])
  
  en.Amplifier(mid=mid, r = r, lab = lab, pos = pos, dtext = dtext) 
}  

## -----------------------------------------------------------------------------

en.Signal <- function (mid,
                    r = 0.03,
                    lab = NULL,
                    pos = 0,
                    dtext = 0.025,
                    ...) {
  user <- par("usr")
  pin <- par("pin")
  sy <- user[4] - user[3]
  sx <- user[2] - user[1]
  ry <- r * sy/sx * pin[1]/pin[2]
  textellipse (mid=mid, radx = r, rady = ry, lab = "")

  if (! is.null(lab))  {
    pp <- tpos(mid, pos, dtext = dtext)
    text(pp$x, pp$y,lab, adj = pp$adj...)
  }   
}

## -----------------------------------------------------------------------------

en.Ground <- function (mid,
                    width = 0.075,
                    length = 0.1,
                    n = 4,
                    dx = 0.2,
                    ...) {
  pp <- mid[2] - length/(n+1)
  segments(mid[1], mid[2], mid[1], pp)

  for (i in 1:n) {
      ww <- width/2 * (1-(i-1) * dx)
      segments(mid[1] - ww, pp, mid[1] + ww, pp)
      pp <- pp - length/(n+1)  
  }
}

