plot.cardiMetacdw <- function(x, y, type = "lattice", scale = TRUE, col.poly = "black", ...) {
  if (is.null(y)) stop("original data must be provided as second argument")

  tt  <- x
  dat <- y

  if (type == "lattice") {
      if (length(tt$weibullfits[[1]]$p) == 4) {
        fweibull <- fweibull4
      } else {
        fweibull <- fweibull6
      }

      xyplot(y ~ x|as.factor(sample), data=dat, ...,
        panel=function(x, y, subscripts, groups, ...) {
          i <- panel.number()
          mres <- tt$weibullfits[[i]]
          smdx <- tt$metares[i, c("tBegin", "tMid", "tEnd")]
          smdy <- fweibull(smdx, mres$p)
          panel.xyplot(x, y)
          xx <- seq(min(x), max(x), length=100)
          yy <- fweibull(xx, mres$p)
          panel.lines(xx, mres$ymax * yy, col="darkgreen", lwd=2)
          #panel.lines(mres$fit$x, mres$fit$f * mres$ymax, col="darkgreen", lwd=2)
          panel.points(smdx, smdy * mres$ymax, col="tomato", pch=16)
       }
      )
  } else {
      if (type == "polygon") {
          par(las=1, font=2, font.lab=2, font.axis=2, mai=c(1.5,1,1,0.5))
          erg <- tt$metares
          count <- nrow(erg)
          if (!(length(col.poly) %in% c(1, count))) stop("length of col.poly does not match number of samples")
          polygoncols <- if(length(col.poly == 1)){rep(col.poly, count)} else {col.poly}
          yrange <- 1:count
          xrange <- floor(seq(min(dat$x), max(dat$x), length=12))
          ymax <- max(erg$yMid)
          ys      <- if(scale == TRUE){erg$yMid / ymax} else {rep(1,count)}
          plot(c(min(xrange),max(xrange)), c(min(yrange), max(yrange) + 1),
            type="n", yaxt="n", ylab="", xlab="Julian Day of Year", ...)
          axis(2, at=yrange, labels=levels(erg$sample))
          abline(h=yrange, lty=3, col=grey(0.7))
          abline(v=xrange, lty=3, col=gray(0.7))
          for (i in 1:count) {
            x <- c(erg$tBegin[i], erg$tMid[i], erg$tEnd[i])
            y <- c(i, i + ys[i], i)
            polygon(x, y, border = NA, col = polygoncols[i])
          }
      } else {stop("plot type unknown")}
  }
}
