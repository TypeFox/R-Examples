plot.spEM <- plot.npEM <- function(x, blocks = NULL, hist=TRUE, addlegend=TRUE,
                                   scale = TRUE, title=NULL, breaks="Sturges", 
                                   ylim=NULL, dens.col, newplot=TRUE, 
                                   pos.legend="topright", cex.legend=1, ...) {
  r <- NCOL(x$data)
  m <- NCOL(x$posteriors)
  blockid <- x$blockid
  if (is.null(blocks)) {
    if(!is.null(blockid)) {
      blocks <- 1:max(blockid)
    } else {
      blocks <- blockid <- 1:r
    }
  }
  ask <- par(ask=(length(unique(blocks))>1))
  ylim.orig <- ylim
  out <- list(x=list(), y=list())
  if (!newplot) {
    hist <- FALSE
  }
  for(i in 1:length(blocks)) {
    coords <- blockid == blocks[i]
    ylim <- ylim.orig
    if (is.null(title)) {
      if (r>1) {
        tt <- paste(which(coords), collapse=",")
        tt <- paste("Coordinate", ifelse(sum(coords)>1, "s ", " "), tt, sep="")
      } else {
        tt <- "Density Curves"
      }
    } else {
      tt <- rep(title,length(blocks))[i]
    }
    dx <- dy <- NULL
    for (j in 1:m) {
      d <- density(x, component=j, block=blocks[i], scale=scale)
      dx <- cbind(dx, d$x)
      dy <- cbind(dy, d$y)
    }
    xx <- as.vector(as.matrix(x$data)[,coords])
    if (is.null(ylim)) {
      ylim=range(dy)
      if (hist) {
        ylim[2] <- max(ylim[2], hist(xx, breaks=breaks, plot=FALSE)$density)
      }
    }
    if (!hist && newplot) {
      pf <- plot # Use plot or hist as plotting fn the 1st time only, then lines
    } else {
      pf <- lines 
    }
    if (hist) {
      hist(xx, breaks=breaks, prob=TRUE, ylim=ylim, main="", ...)
    }
    if (missing(dens.col)) 
      dens.col <- 2:(m+1)
    dens.col <- rep(dens.col, length.out=m)
    for (j in 1:m) {
      pf(dx[,j],dy[,j], type="l", lwd=2, col=dens.col[j], ylim=ylim, ...)
      pf <- lines
    }
    if (addlegend) {
      legend(pos.legend, legend=round(x$lambdahat,3), fill=dens.col, cex=cex.legend)
      out$x[[i]]<-dx
      out$y[[i]]<-dy
    }
    if (newplot) {
      title(main=tt, ...)
    }
  }
  par(ask=ask)
  invisible(out)
}

    

