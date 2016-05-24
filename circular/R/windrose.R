
#############################################################
#                                                           #
#       Original code: Matt Pocernich                       #
#       E-mail: pocernic@rap.ucar.edu                       #
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* # 
# ** Copyright UCAR (c) 1992 - 2004                         #
# ** University Corporation for Atmospheric Research(UCAR)  # 
# ** National Center for Atmospheric Research(NCAR)         #
# ** Research Applications Program(RAP)                     #
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA       #
# ** 2004/28/6 11:31:8                                      #
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* # 
#############################################################

#############################################################
#                                                           #
#   windrose function                                       #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: April, 21, 2011                                   #
#   Version: 0.5                                            #
#                                                           #
#   Copyright (C) 2011 Claudio Agostinelli                  #
#                                                           #
#############################################################

windrose <- function(x, y=NULL, breaks=NULL, bins=12, increment = 10, main='Wind Rose', cir.ind = 0.05, fill.col=NULL, plot.mids=TRUE, mids.size=1.2, osize=0.1, axes=TRUE, ticks=TRUE, tcl=0.025, tcl.text=-0.15, cex=1, digits=2, units=NULL, template=NULL, zero=NULL, rotation=NULL, num.ticks=12, xlim=c(-1.2, 1.2), ylim=c(-1.2, 1.2), uin=NULL, tol=0.04, right=FALSE, shrink=NULL, label.freq=FALSE, calm=c("0", "NA"), ...) {

   calm <- match.arg(calm)
###### internal function used to plot circles
   circles<- function(rad, sector=c(0, 2*pi), lty=2, col="white", border=NA, fill = FALSE) { ## draws circle, radius in units of plot
      ## internal function for windrose
      values <- seq(sector[1], sector[2], by=(sector[2]-sector[1])/360)
      x <- rad*cos(values)
      y <- rad*sin(values)
      if (fill) {
          polygon(x, y, xpd=FALSE, lty=lty, col=col, border=border)
      }
      lines(x, y, col = 1, lty = lty)
   }
####

   if (!is.null(cir.ind) && (cir.ind > 1 | cir.ind <= 0)) {
       cir.ind <- 0.05
       warning("'cir.ind' must be in (0, 1]")
   }  

   if (any(is.null(y))) {
       if (is.data.frame(x) && NCOL(x)==2) {
           y <- x[,2]
           x <- x[,1]
       } else {
           stop("'y' can not be 'NULL' if 'x' is not a dataframe with 2 columns (direction and magnitude) ")
       }
   } else {
       y <- as.vector(y)
       if (length(y)!=length(x)) stop("'x' and 'y' must have the same length")
       if (is.data.frame(x)) {
           if (NCOL(x) > 1) x <- x[,1]
       } 
   }

### rm calms ###
#### count calms
#### following NWS protocol, when dir == 0 weather calm, we allows also to set cam==NA
   
   no <- length(x)
   if (calm=="NA") {
       calmcalm <- sum(is.na(x))
       notcalm <- !is.na(x)
   } else {
       calmcalm <- sum(x == calm)
       notcalm <- x!=calm
   }
   x <- x[notcalm]
   y <- y[notcalm]
  
   # Handling missing values in any case
   ok <- complete.cases(x, y)
   x <- x[ok]
   y <- y[ok]
    
   if (length(y)==0) {
       warning("No observations (at least after removing missing values and calm winds)")
       return(NULL)
   }
   
   result <- list()
   result$x <- x
   result$y <- y
   xcircularp <- attr(as.circular(x), "circularp")
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
   } else {
      if (is.null(zero))
         zero <- xcircularp$zero
      if (is.null(rotation))
         rotation <- xcircularp$rotation
   }

   op <- par(mar = c(1,1,2,1))
   mai <- par("mai") 
   on.exit(par(op))

   midx <- 0.5 * (xlim[2] + xlim[1])
   xlim <- midx + (1 + tol) * 0.5 * c(-1, 1) * (xlim[2] - xlim[1])
   midy <- 0.5 * (ylim[2] + ylim[1])
   ylim <- midy + (1 + tol) * 0.5 * c(-1, 1) * (ylim[2] - ylim[1])
   oldpin <- par("pin") - c(mai[2]+mai[4], mai[1]+mai[3])
   xuin <- oxuin <- oldpin[1]/diff(xlim)
   yuin <- oyuin <- oldpin[2]/diff(ylim)
   if (is.null(uin)) {
       if (yuin > xuin) xuin <- yuin
       else yuin <- xuin
   } else {
       if (length(uin) == 1) uin <- uin * c(1, 1)
       if (any(c(xuin, yuin) < uin)) stop("uin is too large to fit plot in")
       xuin <- uin[1]; yuin <- uin[2]
   }    
   xlim <- midx + oxuin/xuin * c(-1, 1) * diff(xlim) * 0.5
   ylim <- midy + oyuin/yuin * c(-1, 1) * diff(ylim) * 0.5
       
   if (any(is.null(breaks))) {
     step <- 2*pi/bins
     breaks <- circular(seq(0, 2*pi, by=step), units="radians")
   } else {
       breaks <- as.circular(breaks)
   }
   breaks <- conversion.circular(breaks, units="radians", zero=0, rotation="counter", modulo="2pi")
   attr(breaks, "class") <- attr(breaks, "circularp") <-  NULL
   if (template=="clock12") { ### added for clock12
     breaks <- 2*breaks
     breaks <- breaks%%(2*pi)
   }
   breaks <- sort(unique(breaks))
   if (breaks[1]!=0) {
       breaks <- c(breaks[length(breaks)]-2*pi, breaks)
   } else {
       breaks <- c(breaks, 2*pi)
   }   
   bins <- length(breaks)-1
   step <- diff(breaks) # the step for the breaks which include zero degrees is the first one
   
   x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
   attr(x, "class") <- attr(x, "circularp") <-  NULL
   
   if (template=="clock12") { ### added for clock12
     x <- 2*x
     x <- x%%(2*pi)
   }
   x[x >= breaks[bins+1]] <- x[x >= breaks[bins+1]]-2*pi
 
###   x[(x+step/2)%%(2*pi)==0] <- 2*pi # Values like 359 go to Sector 0
   plot(c(-1.2,1.2), c(-1.2,1.2), xlab='', ylab='', main=main, xaxt='n', yaxt='n', pch=' ', xlim=xlim, ylim=ylim)
   counts <- hist.default(x, breaks=breaks ,plot=FALSE, right=right)$counts #use hist for
   mids <- breaks[1:bins]+step/2 # midpoints
   if (plot.mids) {
       for (i in 1:bins) { 
            lines(c(0, mids.size*cos(mids[i])), c(0, mids.size*sin(mids[i])), lty=2) 
       }
   }

#######  lines
   maxlength.orig <- sqrt(max(counts/step)/length(x)+osize^2)
   if (is.null(shrink)) {
       maxlength <- shrink <- maxlength.orig
   } else {
       maxlength <- shrink
   }
   J <- ceiling(max(y, na.rm = TRUE)/increment)  ## should deal with NA elsewhere
   if (any(is.null(fill.col))) {
       if (J%%2 == 0) {
           fill.col <- c("blue", "red")
       } else {
           fill.col <-  c("red", "blue")
       }
   }
   fill.col <- rep(fill.col, length.out=bins)

   OUT <- matrix(NA, nrow = J, ncol = bins)
   OUT[J,] <- counts

   for(j in J:1){  
      data1<- x[y <= j*increment]
      OUT[j,] <- counts <- hist.default(data1, breaks=breaks, plot=FALSE, right=right)$counts #use hist for
      for (i in 1:bins) {
           w1 <- breaks[i]  ## in radians, the locations of the upper and lower lines
           w2 <- breaks[i+1]

           if (counts[i]) {
               rad <- sqrt(counts[i]/(step[i]*length(x)) + osize^2)/maxlength
           } else {
               rad <- 0
           }
           xx <- rad*c(0,cos(w1),cos(w2),0) ## increase length by percent equal to bin with 
           yy <- rad*c(0,sin(w1),sin(w2),0)
           polygon(xx, yy, xpd=FALSE, col = fill.col[j], border=NA)
   
           lines(xx[1:2], yy[1:2])
           lines(xx[3:4], yy[3:4])
           circles(rad=rad, sector=c(w1, w2), fill=TRUE, lty=1, col=fill.col[j], border=NA) 
      } ## close for i
   } ## close J loop

   m <- dim(OUT)[1]
   new <- OUT
   if (m > 1) {
       for (i in 2:m) { 
            new[i,]<- OUT[i,] - OUT[i-1,]
       }
   }

   ##  circles
   if (!is.null(cir.ind)) {
       equalstep <- max(abs(diff(step))) <= 10*.Machine$double.eps
       max.plt <- maxlength.orig^2 - osize^2
       if (equalstep & label.freq) max.plt <- max.plt*step[1]
       cir.ind <- min(cir.ind, max.plt)
       max.plt <- floor(max.plt/cir.ind)*cir.ind  ## sets max plotted area      
       max.plt <- seq(cir.ind, max.plt, by = cir.ind)
       if (equalstep & label.freq) {
           rad <- sqrt(max.plt/step[1]+osize^2)/maxlength
       } else {
           rad <- sqrt(max.plt+osize^2)/maxlength     
       }
       if (equalstep & label.freq) {
           text(0, rad, paste(round(max.plt * 100, digits=digits), "%", sep = ""), cex = 0.9*cex, pos = 3, font = 3, offset=0.2)
       } else {
           text(0, rad, paste(round(max.plt, digits=digits), sep = ""), cex = 0.9*cex, pos = 3, font = 3, offset=0.2)
       }
       for (i in 1:length(rad)) {
            circles(rad[i])
       } # close circle
   } else {
       circles(1)
   }

   circles(osize, fill = TRUE)
              
   if (axes) {
     axis.circular(at=NULL, labels=NULL, units=units, template=template, modulo="2pi", zero=zero, rotation=rotation, tick=ticks, cex=cex, tcl=tcl, tcl.text=tcl.text, digits=digits)
   }

   if (axes==FALSE & ticks) {
       at <- (0:num.ticks)/num.ticks*2*pi
       if (rotation=="clock") at <- -at
       at <- at + zero
       ticks.circular(circular(x=at, type="angles", units="radians", modulo="asis", zero=zero, rotation=rotation), tcl=tcl)
   }
   
   OUT <- round(new/sum(OUT[m,]),3)
   colnamesout <- rep("", bins)
   breaks <- conversion.circular(circular(breaks), units=units)
   mids <- conversion.circular(circular(mids), units=units)
   for (i in 1:bins) {
        if (right) {
            colnamesout[i] <- paste("(", round(breaks[i], digits=digits), ", ", round(breaks[i+1], digits=digits), "]", sep="")
        } else {
            colnamesout[i] <- paste("[", round(breaks[i], digits=digits), ", ", round(breaks[i+1], digits=digits), ")", sep="")
        }
   }
   colnames(OUT) <- colnamesout
   rownames(OUT) <- paste( "(", 0:(J-1)*increment,  ",", 1:J * increment, "]", sep = "")

   result$table <- OUT
   result$number.obs <- no
   result$number.calm <- calmcalm
   result$breaks <- breaks
   result$mids <- mids
   result$shrink <- shrink
   result$call <- match.call()

   invisible(result)
}
