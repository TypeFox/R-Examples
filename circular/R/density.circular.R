#############################################################
#                                                           #
#   density.circular function                               #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   date: February, 14, 2013                                #
#   Copyright (C) 2013 Claudio Agostinelli                  #
#                                                           #
#   Version 0.3-1                                           #
#                                                           #
#############################################################

density.circular <- function(x, z=NULL, bw, adjust = 1, type = c("K", "L"), kernel= c("vonmises", "wrappednormal"), na.rm = FALSE, from=circular(0), to=circular(2*pi), n=512, K=NULL, min.k=10, control.circular=list(), ...) {

    name <- deparse(substitute(x))
    data <- x

    if (!is.numeric(from))
        stop("argument 'from' must be numeric")      
    if (!is.numeric(to))
        stop("argument 'to' must be numeric")      
    if (!is.finite(from)) 
        stop("non-finite `from'")
    if (!is.finite(to)) 
        stop("non-finite `to'")
    if (!is.numeric(n))
        stop("argument 'n' must be numeric")
    n <- round(n)
    if (n <=0)
         stop("argument 'n' must be integer and positive")         
    if (!is.numeric(x)) 
        stop("argument 'x' must be numeric")
    if (!is.null(z) && is.circular(z)) {
       datacircularp <- circularp(z)
    } else if (is.circular(x))
              datacircularp <- circularp(x)
    else {
       datacircularp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
    }
    dc <- control.circular
    if (is.null(dc$type))
       dc$type <- datacircularp$type
    if (is.null(dc$units))
       dc$units <- datacircularp$units
    if (is.null(dc$template))
       dc$template <- datacircularp$template
    if (is.null(dc$modulo))
       dc$modulo <- datacircularp$modulo
    if (is.null(dc$zero))
       dc$zero <- datacircularp$zero
    if (is.null(dc$rotation))
       dc$rotation <- datacircularp$rotation
    if (dc$modulo=="pi")
      stop("The function does not work yet for modulo='pi'")

    x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
    attr(x, "class") <- attr(x, "circularp") <- NULL
    from <- conversion.circular(from, units="radians", zero=0, rotation="counter")
    attr(from, "class") <- attr(from, "circularp") <- NULL
    to <- conversion.circular(to, units="radians", zero=0, rotation="counter")
    attr(to, "class") <- attr(to, "circularp") <- NULL

    kernel <- match.arg(kernel)

    x <- as.vector(x)
    x.na <- is.na(x)
    if (any(x.na)) {
        if (na.rm) 
            x <- x[!x.na]
        else stop("x contains missing values")
    }
    
    x.finite <- is.finite(x)
    if (any(!x.finite)) {
        x <- x[x.finite]
    }
    nx <- length(x)
    
    if (is.null(z)) {
        z <- circular(seq(from=from, to=to, length=n))
    } else {
        if (!is.numeric(z))
           stop("argument 'z' must be numeric")
        namez <- deparse(substitute(z))
        z.na <- is.na(z)
        if (any(z.na)) {
            if (na.rm) {
                z <- z[!z.na]
            } else {
                stop("z contains missing values")
            }
        }
        z.finite <- is.finite(z)
        if (any(!z.finite)) {
            z <- z[z.finite]
        }
    }
    zz <- conversion.circular(z, dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
    z <- conversion.circular(z, units="radians", zero=0, rotation="counter")
    attr(z, "class") <- attr(z, "circularp") <- NULL
    z <- as.vector(z)
    
    bw <- adjust * bw
    if (!is.numeric(bw))
        stop("argument 'bw' and 'adjust' must be numeric")        
    if (!is.finite(bw)) 
        stop("non-finite `bw'")
    if (bw <= 0) 
        stop("`bw' is not positive.")
    
    y <- DensityCircularRad(x=x, z=z, bw=bw, kernel=kernel, K=K, min.k=min.k)

    structure(list(data = data, x = zz, y = y, bw = bw, n = nx, kernel=kernel, call = match.call(), data.name=name, has.na = FALSE), class = "density.circular")
} 

DensityCircularRad <- function(x, z, bw, kernel, K=NULL, min.k=10) {
   nx <- length(x)
   if (kernel=="vonmises") {
       y <- sapply(z, DvonmisesRad, mu=x, kappa=bw, log=FALSE)
   } else if (kernel=="wrappednormal") {
       rho <- exp(-bw^2/2)
       y <- sapply(z, DwrappednormalRad, mu=x, rho=rho, K=K, min.k=min.k)
   } else {
       stop("other kernels not implemented yet")
   }
   y <- apply(y, 2, sum)/nx
   return(y)
}

#############################################################
#                                                           #
#   plot.density.circular function                          #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 22, 2011                                     #
#   Copyright (C) 2011 Claudio Agostinelli                  #
#                                                           #
#   Version 0.5-12                                           #
#                                                           #
#############################################################

plot.density.circular <- function(x, main = NULL, sub=NULL, xlab = NULL, ylab ="Density circular", type = "l", zero.line = TRUE, points.plot=FALSE, points.col=1, points.pch=1, points.cex=1, plot.type = c("circle", "line"), axes=TRUE, ticks=FALSE, bins=NULL, offset=1, shrink=1, tcl=0.025,  tcl.text=0.125, sep=0.025, tol = 0.04, digits=2, cex=1, uin=NULL, xlim=NULL, ylim=NULL, join=FALSE, nosort=FALSE, units=NULL, template=NULL, zero=NULL, rotation=NULL, control.circle=circle.control(), ...) {
   xcircularp <- attr(x$x, "circularp")
   if (is.null(xcircularp))
      stop("the component 'x' of the object must be of class circular")
##   type <- xcircularp$type
   modulo <- xcircularp$modulo
   if (is.null(units)) 
      units <- xcircularp$units
   if (is.null(template))
      template <- xcircularp$template
   if (template=="geographics" | template=="clock24") {
      zero <- pi/2
      rotation <- "clock"
   } else if (template=="clock12") {
      zero <- pi/2
      rotation <- "clock"
      modulo <- "pi"
   } else {
      if (is.null(zero))
         zero <- xcircularp$zero
      if (is.null(rotation))
         rotation <- xcircularp$rotation
   }
   next.points <- 0
   if (template=="clock12") {
     x$x <- 2*x$x
     x$data <- 2*x$data
   } 
   x$x <- conversion.circular(x$x, units="radians", modulo="2pi")
   x$data <- conversion.circular(x$data, units="radians", modulo="2pi")
   attr(x$x, "class") <- attr(x$x, "circularp") <- NULL
   attr(x$data, "class") <- attr(x$data, "circularp") <- NULL

   plot.type <- match.arg(plot.type)

   if (is.null(xlab))
      xlab <- paste("N =", x$n, "  Bandwidth =", formatC(x$bw), " Unit =", units)
   if (is.null(main))
      main <- deparse(x$call)
#### as scatter plot
   if (plot.type == "line") {
      if (units=='degrees') {
        x$x <- x$x/pi*180
        x$data <- x$data/pi*180
      }
      if (units=='hours') {
        x$x <- x$x/pi*12
        x$data <- x$data/pi*12
      }
      if (is.null(xlim))
         xlim <- range(c(x$x, x$data))
      if (is.null(ylim)) {
         ylim <- range(x$y)
         if (points.plot)
           ylim[1] <- ylim[1]-0.04*points.cex
      }
      
      xorder <- order(x$x)
      x$x <- x$x[xorder]
      x$y <- x$y[xorder]
      
      plot.default(x, type = type, xlim=xlim, ylim=ylim, main=main, xlab=xlab, ylab=ylab, axes=axes, ...)
      if (zero.line) 
         abline(h = 0, lwd = 0.2, col = "gray")
      if (points.plot)
         points(x$data, rep(ylim[1]+0.02*points.cex, length(x$data)), col=points.col, pch=points.pch, cex=points.cex)
      return(NULL)
   } else {
#### as circular plot
      if (is.null(xlim))
         xlim <- c(-1, 1)
      if (is.null(ylim))
         ylim <- c(-1, 1)
      if (is.null(bins)) {
         bins <- NROW(x)
      } else {
         bins <- round(bins)
         if (bins<=0)
            stop("'bins' must be non negative")
      }
    
      CirclePlotRad(xlim=xlim, ylim=ylim, uin=uin, shrink=shrink, tol=tol, main=main, sub=sub, xlab=xlab, ylab=ylab, control.circle=control.circle)

      if (!is.logical(ticks))
        stop("ticks must be logical")
      
      if (axes) {
	 axis.circular(at=NULL, labels=NULL, units=units, template=template, modulo=modulo, zero=zero, rotation=rotation, tick=ticks, cex=cex, tcl=tcl, tcl.text=tcl.text, digits=digits)
      }
      if (axes==FALSE & ticks) {
         at <- circular((0:bins)/bins*2*pi, zero=zero, rotation=rotation)
         ticks.circular(at, tcl=tcl)
      }

      if (rotation=="clock")
         x$x <- -x$x
      x$x <- x$x+zero
      x$x <- x$x%%(2*pi)
      ll <- LinesCircularRad(x=x$x, y=x$y, join=join, nosort=nosort, offset=offset, shrink=shrink, ...)

      if (points.plot) {
         next.points <- sep
         if (rotation=="clock")
            x$data <- -x$data
         x$data <- x$data+zero
         x$data <- x$data%%(2*pi)
         PointsCircularRad(x$data, bins, FALSE, points.col, points.pch, 1, 1, sep, next.points, shrink, points.cex)
      }
      return(invisible(list(x=ll$x, y=ll$y, zero=zero, rotation=rotation, next.points=next.points)))
   }
}


#############################################################
#                                                           #
#   lines.density.circular function                         #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 5, 2011                                     #
#   Copyright (C) 2010 Claudio Agostinelli                  #
#                                                           #
#   Version 0.5-1                                           #
#                                                           #
#############################################################

lines.density.circular <- function(x, type = "l", zero.line = TRUE, points.plot=FALSE, points.col=1, points.pch=1, points.cex=1, plot.type = c("circle", "line"), bins=NULL, offset=1, shrink=1, tcl=0.025, sep=0.025, join=TRUE, nosort=FALSE, plot.info=NULL, zero=NULL, rotation=NULL, ...) {
   xcircularp <- attr(x$x, "circularp")
   if (is.null(xcircularp))
      stop("the component 'x' of the object must be of class circular")
##   type <- xcircularp$type
   modulo <- xcircularp$modulo
   units <- xcircularp$units
   template <- xcircularp$template
   if (is.null(plot.info)) {
      if (is.null(zero)) {
        if (template=="geographics" | template=="clock24" | template=="clock12")
           zero <- pi/2
        else
           zero <- xcircularp$zero
      }
      if (is.null(rotation)) {
        if (template=="geographics" | template=="clock24" | template=="clock12")
          rotation <- "clock"
        else
          rotation <- xcircularp$rotation
      }
      next.points <- 0
   } else {
      zero <- plot.info$zero
      rotation <- plot.info$rotation
      next.points <- plot.info$next.points
   }
   
   if (template=="clock12") {
     x$x <- 2*x$x
     x$data <- 2*x$data
   } 
   
   x$x <- conversion.circular(x$x, units="radians")
   x$data <- conversion.circular(x$data, units="radians")
   attr(x$x, "circularp") <- attr(x$x, "class") <- NULL
   attr(x$data, "circularp") <- attr(x$data, "class") <- NULL

   ll <- list()
   
   plot.type <- match.arg(plot.type)
   if (is.null(bins)) {
       bins <- NROW(x)
   } else {
       bins <- round(bins)
       if (bins<=0)
          stop("bins must be non negative")
   }
    
    if (plot.type == "line") {
        if (units=='degrees') {
          x$x <- x$x/pi*180
          x$data <- x$data/pi*180
        }
        if (units=='hours') {
          x$x <- x$x/pi*12
          x$data <- x$data/pi*12
        }
        xorder <- order(x$x)
        x$x <- x$x[xorder]
        x$y <- x$y[xorder] 
        lines.default(x, type = type, ...)
        if (zero.line) 
            abline(h = 0, lwd = 0.2, col = "gray")
        if (points.plot)
            points.default(x$data, rep(min(x$y)-0.02*points.cex, length(x$data)), col=points.col, pch=points.pch)
    } else {
      if (rotation=="clock")
         x$x <- -x$x
      x$x <- x$x+zero
      x$x <- x$x%%(2*pi)
      ll <- LinesCircularRad(x=x$x, y=x$y, join=join, nosort=nosort, offset=offset, shrink=shrink, ...)
      if (points.plot) {
         if (rotation=="clock")
            x$data <- -x$data
         x$data <- x$data+zero
         x$data <- x$data%%(2*pi)
         next.points <- next.points+sep
         PointsCircularRad(x$data, bins, FALSE, points.col, points.pch, 1, 1, sep, next.points, shrink, points.cex)
      }
    }
      
    return(invisible(list(x=ll$x, y=ll$y, zero=zero, rotation=rotation, next.points=next.points)))
}


#############################################################
#                                                           #
#   print.density.circular function                         #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 23, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#                                                           #
#############################################################

print.density.circular <- function(x, digits=NULL, ...)
{
    cat("\nCall:\n\t",deparse(x$call),
	"\n\nData: ",x$data.name," (",x$n," obs.);",
	"\tBandwidth 'bw' = ",formatC(x$bw,digits=digits), "\n\n",sep="")
    print(summary(as.data.frame(x[c("x","y")])), digits=digits, ...)
    invisible(x)
}

