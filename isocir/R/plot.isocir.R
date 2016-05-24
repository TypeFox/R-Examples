
"plot.isocir" <- function(x, option=c("CIRE","cirmeans"), cex=1, stack=TRUE, axes=TRUE, sep=0.025, shrink=1, bins=300, ticks=FALSE, tcl=0.025, tcl.text=0.125, col=NULL, tol=0.04, uin=NULL, xlim=c(-1, 1), ylim=c(-1, 1), digits=2, units=NULL, template=NULL, zero=NULL, rotation=NULL, main=NULL, sub=NULL, xlab="", ylab="", control.circle=circle.control(), ...) {

# DATE WRITTEN: 1 Dic 2010          LAST REVISED:  09 Jan 2012
# AUTHOR: Sandra Barragan based on the code of Ulric Lund and Claudio Agostinelli
# DESCRIPTION: This is an auxiliary function to enable class isocir to be plotted.
# REFERENCE: It is mainly based on some functions of the package circular (by Ulric Lund and Claudio Agostinelli) but with some changes.
# SEE ALSO: CIRE, cond.test.


CirclePlotRad <- function(xlim=c(-1,1), ylim=c(-1,1), uin=NULL, shrink=1, tol=0.04, main=NULL, sub=NULL, xlab=NULL, ylab=NULL, control.circle=circle.control()) {
   xlim <- shrink * xlim
   ylim <- shrink * ylim
   midx <- 0.5 * (xlim[2] + xlim[1])
   xlim <- midx + (1 + tol) * 0.5 * c(-1, 1) * (xlim[2] - xlim[1])
   midy <- 0.5 * (ylim[2] + ylim[1])
   ylim <- midy + (1 + tol) * 0.5 * c(-1, 1) * (ylim[2] - ylim[1])
   oldpin <- par("pin")
   xuin <- oxuin <- oldpin[1]/diff(xlim)
   yuin <- oyuin <- oldpin[2]/diff(ylim)
   if (is.null(uin)) {
      if (yuin > xuin)
         yuin <- xuin
      else
         xuin <- yuin
   } else {
      if (length(uin) == 1)
         uin <- uin * c(1, 1)
      if (any(c(xuin, yuin) < uin))
         stop("uin is too large to fit plot in")
      xuin <- uin[1]; yuin <- uin[2]
   }
   xlim <- midx + oxuin/xuin * c(-1, 1) * diff(xlim) * 0.5
   ylim <- midy + oyuin/yuin * c(-1, 1) * diff(ylim) * 0.5
   n <- control.circle$n
   x <- cos(seq(0, 2 * pi, length = n))
   y <- sin(seq(0, 2 * pi, length = n))
   axes <- FALSE
   log <- ""
   xaxs <- "i"
   yaxs <- "i"
   ann <-  par("ann")
   frame.plot <- axes
   panel.first <- NULL 
   panel.last <- NULL
   asp <- NA
   plot.default(x=x, y=y, type=control.circle$type, xlim=xlim, ylim=ylim, log="", main=main, sub=sub, xlab=xlab, ylab=ylab, ann=ann, axes=axes, frame.plot=frame.plot, panel.first=panel.first, panel.last=panel.last, asp=asp, col=control.circle$col, bg=control.circle$bg, pch=control.circle$pch, cex=control.circle$cex, lty=control.circle$lty, lwd=control.circle$lwd)
}

circle.control <- function(n=1000, type='l', col=1, bg=par('bg'), pch=1, cex=1, lty=1, lwd=1) {
  x <- list(n=n, type=type, col=col, bg=bg, pch=pch, cex=cex, lty=lty, lwd=lwd)
  return(x)
}

PointsCircularmRad <- function(x, bins, stack, col, pch, iseries, nseries, sep, next.points, shrink, cex, ...) {
#### x musts be in modulo 2pi  
   if (!stack) {
      z <- cos(x)
      y <- sin(x)
      r <- 1+((iseries-1)*sep+next.points)*shrink
      points.default(z*r, y*r, cex=cex, pch=pch, col = col[iseries], ...)
   } else {
      x[x >= 2*pi] <- 2*pi-4*.Machine$double.eps
      arc <- (2 * pi)/bins
      pos.bins <- ((1:nseries)-1/2)*arc/nseries-arc/2
      breaks <- seq(0,2*pi,length.out=(bins+1))
      bins.count <- hist.default(x, breaks=breaks, plot=FALSE, right=TRUE)$counts
      mids <- seq(arc/2, 2 * pi - pi/bins, length = bins) + pos.bins[iseries]
      index <- cex*sep
	l<-0
      for (i in 1:bins) {
         if (bins.count[i] != 0) {
            for (j in 0:(bins.count[i] - 1)) {
               r <- 1 + j * index
               z <- r * cos(mids[i])
               y <- r * sin(mids[i])
			l<-l+1
               points.default(z, y, cex=cex, pch=pch[l], col=col[iseries], ...)
            } # for j
         } # if
      } # for i
   } # else
} # function Points


   if(class(x)!= "isocir"){stop("The argument x must be a object of class isocir")}
   option <- match.arg(option)
   if(option=="CIRE"){
    x <- unlist(x$CIRE)
    main="Circular Isotonic Regression Estimator"
    col=2
   }
   else if(option=="cirmeans"){
    x <- unlist(x$cirmeans)
    main="Unrestricted Estimator"
    col=1
   }
   else{stop("option must be CIRE or cirmeans")}

   pch <- c(1:length(x))

   if (is.matrix(x) | is.data.frame(x)) {
      nseries <- ncol(x)
   } else {
      nseries <- 1
   }
   xx <- as.data.frame(x)
  
   xcircularp <- attr(suppressWarnings(as.circular(xx[,1])), "circularp")
   type <- xcircularp$type
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
   
   CirclePlotRad(xlim=xlim, ylim=ylim, uin=uin, shrink=shrink, tol=tol, main=main, sub=sub, xlab=xlab, ylab=ylab, control.circle=control.circle)
    
   if (is.null(bins)) {
      bins <- NROW(x)
   } else {
      bins <- round(bins)
      if (bins<=0)
         stop("bins must be non negative")
   }

   if (is.null(col)) {
      col <- seq(nseries)
   } else {
      if (length(col)!=nseries) {
         col <- rep(col, nseries)[1:nseries]
      }
   }
if(length(pch)==1){pch <- rep(pch, nseries, length.out=nseries)}

   if (axes) {
      axis.circular(units=units, template=template, zero=zero, rotation=rotation, digits=digits, cex=cex, tcl=tcl, tcl.text=tcl.text)
   }
   
   if (!is.logical(ticks))
      stop("ticks must be logical")
   
   if (ticks) {
      at <- circular((0:bins)/bins*2*pi, zero=zero, rotation=rotation)
      ticks.circular(at, tcl=tcl)
   }

   for (iseries in 1:nseries) {      
      x <- xx[,iseries]
      x <- na.omit(x)
      n <- length(x)      
      if (n) {
         x <- suppressWarnings(conversion.circular(x, units="radians", modulo=modulo))
         attr(x, "circularp") <- attr(x, "class") <- NULL
         if (rotation=="clock")
            x <- -x
         x <- x+zero
         if (template=="clock12")
           x <- 2*x
         x <- x%%(2*pi)
         PointsCircularmRad(x, bins, stack, col, pch, iseries, nseries, sep, 0, shrink, cex, ...)
      }
   }
return(invisible(list(zero=zero, rotation=rotation, next.points=nseries*sep)))
}