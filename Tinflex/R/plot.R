#############################################################################
##
##  Plot hat and squeeze.
##
#############################################################################


plot.Tinflex <- function(x, from, to, is.trans=FALSE, ...) {
  ## ------------------------------------------------------------------------
  ## S3 method for plotting class 'Tinflex'.
  ## ------------------------------------------------------------------------
  ##   x        ... S3 object of class 'Tinflex'
  ##   from, to ... range over which the function will be plotted
  ##   is.trans ... whether transformed scale is used
  ##   ...      ... graphical parameters
  ## ------------------------------------------------------------------------

  ## Check arguments.
  if (missing(from) || missing(to))
    stop ("argument 'from' or 'to' is missing, with no default")
  
  ## Get parameters for hat and squeeze.
  ivs <- x$ivs
  n.ivs <- ncol(ivs)-1
  
  ## Plot (transformed) pdf.
  if (isTRUE(is.trans)) {
    Tpdf <- function(z) { Tf(x$lpdf, ivs["c",1], z) }
    plot(Tpdf, from=from, to=to, n=501, col="blue", ...)
  }
  else {
    pdf <- function(z) { exp(x$lpdf(z)) }
    plot(pdf, from=from, to=to, n=501, col="blue", ...)
  }

  ## Plot hat and squeeze.
  ## The graphs are plotted piecewise for each of the intervals.
  for (i in 1:n.ivs) {

    ## Create x values.
    xb <- c(ivs["x",i],ivs["x",i+1])
    xb <- pmax(xb, from)
    xb <- pmin(xb, to)
    
    n <- max(10, floor(500*(xb[2]-xb[1])/(to-from)))
    h <- (0:n)/n
    x <- (1-h)*xb[1]+h*xb[2]

    ## Plot (piece of) squeeze.
    squeeze <- function(z) { ivs["sq.a",i] + ivs["sq.b",i]*(z-ivs["sq.y",i])}
    if (! (is.na(ivs["sq.a",i]) || is.na(ivs["sq.a",i]))) {
      if (isTRUE(is.trans)) {
        y <- squeeze(xb)
        lines(xb,y, col="darkgreen", lwd="2")
      }
      else {
        y <- Tinv(ivs["c",i], squeeze(x))
        lines(x,y, col="darkgreen", lwd="2")
      }
    }

    ## Plot (piece of) hat.
    hat <- function(z) { ivs["ht.a",i] + ivs["ht.b",i]*(z-ivs["ht.y",i])}
    if (isTRUE(is.trans)) {
      y <- hat(xb)
      lines(xb,y, col="red", lwd="2")
    }
    else {
      y <- Tinv(ivs["c",i], hat(x))
      lines(x,y, col="red", lwd="2")
    }

  }
}

## --------------------------------------------------------------------------
