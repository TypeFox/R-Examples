"plotradarlmr" <-
function(lmom, num.axis=4, plot=TRUE, points=FALSE, poly=TRUE, tag=NA,
               title="L-moment Ratio Radar Plot", make.zero.axis=FALSE,
               minrat=NULL, maxrat=NULL, theomins=TRUE, rot=0,
               labadj=1.2, lengthadj=1.75, offsetadj=0.25, scaleadj=2.2,
               axis.control  = list(col=1, lty=2, lwd=0.5, axis.cex=0.75, lab.cex=0.95),
               point.control = list(col=8, lwd=0.5, pch=16),
               poly.control  = list(col=rgb(0,0,0,.1), border=1, lty=1, lwd=1), ...) {

   if(! is.null(lmom)) {
      if((length(lmom$ratios) - 2) < num.axis) { 
         num.axis <- (length(lmom$ratios) - 2)
      }
      rats <- lmom$ratios[3:(num.axis+2)]
      if(length(rats) < 3) {
         warning("The order of L-moment ratios must be greater than 5 (at least through L5)")
         return()
      }
   } else if(make.zero.axis) {
      rats <- rep(0, num.axis)
   } else {
      rats <- rep(NA,num.axis)
   }

   if(is.null(minrat)) {
      minrat <- minrat.text <- rep(-1, length(rats))
   } else if(length(minrat) == 1) {
      minrat      <- rep(as.numeric(minrat),   length(rats))
      minrat.text <- rep(as.character(minrat), length(rats))
   } else {
      if(length(minrat) != num.axis) {
         warning("Inconsistency of user-level minrat to num.axis provided")
         return()
      }
      minrat.text <- minrat
   }

   if(is.null(maxrat)) {
      maxrat <- maxrat.text <- rep( 1, length(rats))
   } else if(length(maxrat) == 1) {
      maxrat      <- rep(as.numeric(maxrat),   length(rats))
      maxrat.text <- rep(as.character(maxrat), length(rats))
   } else {
      if(length(maxrat) != num.axis) {
         warning("Inconsistency of user-level maxrat to num.axis provided")
         return()
      }
      maxrat.text <- maxrat
   }

   if(theomins) {
      minrat.text[2] <- "-1/4"; minrat[2] <- -1/4
      if(num.axis > 3) {
         minrat.text[4] <- "-1/6"; minrat[4] <- -1/6
      }
   }
   rotation  <- rot

   theta    <- 360/num.axis
   thetas   <- seq(0, 360, by=theta) + rotation
   thetas.r <- 2*pi*thetas/360
   if(length(thetas.r) > num.axis) {
      thetas.r <- thetas.r[1:(length(thetas.r)-1)]
   } else if(length(thetas.r) < num.axis) {
      warning("Computed angles of the axes have incorrect length to provided num.axis")
      return()
   }

   max.m.min <- (maxrat - minrat)/lengthadj

   axes <- FALSE
   false.plot.stub <- labadj*c(-scaleadj*labadj*lengthadj, scaleadj*labadj*lengthadj)
   if(plot) {
       plot(false.plot.stub,false.plot.stub, type="n", xaxt="n", yaxt="n",
            bty="n", xlab="", ylab="")
       axes <- TRUE
   }
   x <- y <- xo1 <- yo1 <- xo2 <- yo2 <- xz <- yz <- vector(mode="numeric")
   for(i in 1:length(thetas.r)) {
      phi <- thetas.r[i]; cosphi <- cos(phi); sinphi <- sin(phi)
      if(axes) {
         leg1 <- abs(minrat[i]) + offsetadj
         legz <- (       0  - minrat[i]) / max.m.min[i] + leg1
         leg2 <- (maxrat[i] - minrat[i]) / max.m.min[i] + leg1
         #print(c(leg1, legz, leg2))
         leg1 <- lengthadj*leg1; leg2 <- lengthadj*leg2; legz <- lengthadj*legz
         xo1[i] <- leg1*cosphi; yo1[i] <- leg1*sinphi
         xo2[i] <- leg2*cosphi; yo2[i] <- leg2*sinphi
          xz[i] <- legz*cosphi;  yz[i] <- legz*sinphi
         arrows(xo1[i],yo1[i],xo2[i],yo2[i], length = 0, angle = theta/8,
                  col=axis.control$col, lty=1, lwd=axis.control$lwd)
         if(axis.control$axis.cex != 0) {
            text(xo1[i],yo1[i], minrat.text[i], cex=axis.control$axis.cex)
            text(xo2[i],yo2[i], maxrat[i], cex=axis.control$axis.cex)
         }
         if(axis.control$lab.cex != 0) {
            text(labadj*leg2*cosphi, labadj*leg2*sinphi, paste("Tau",(i+2), sep=""),
                 cex=axis.control$lab.cex)
         }
      }
      leg <- (rats[i] - minrat[i]) / max.m.min[i] + abs(minrat[i]) + offsetadj
      x[i] <- lengthadj*leg*cosphi; y[i] <- lengthadj*leg*sinphi
   }
   if(axes) {
      xo1.cycle <- c(xo1,xo1[1]); yo1.cycle <- c(yo1,yo1[1]) # the minimum origin
      xo2.cycle <- c(xo2,xo2[1]); yo2.cycle <- c(yo2,yo2[1]) # the maximum origin
       xz.cycle <- c(xz,  xz[1]);  yz.cycle <- c(yz,  yz[1]) # the zero origin

      polygon(xo1.cycle, yo1.cycle, border=axis.control$col, lty=1, lwd=axis.control$lwd)
      polygon(xo2.cycle, yo2.cycle, border=axis.control$col, lty=1, lwd=axis.control$lwd)
      polygon( xz.cycle,  yz.cycle, border=axis.control$col, lty=axis.control$lty, lwd=axis.control$lwd)
   }
   x.cycle <- c(x, x[1]); y.cycle <- c(y, y[1])
   if(poly) {
      polygon(x.cycle, y.cycle,
              col=poly.control$col, border=poly.control$border, lty=poly.control$lty, lwd=poly.control$lwd)
   }
   if(points) {
      points(x,y,
             col=point.control$col, pch=point.control$pch, lwd=point.control$lwd)
   }
   if(plot) {
      if(! is.na(title)) text(0, false.plot.stub[2], title, ...)
      if(! is.na(tag))   text(0, 0, tag, ...)
   }
}

