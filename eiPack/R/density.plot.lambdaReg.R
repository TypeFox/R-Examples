density.plot.lambdaReg <- function(x, by = "column", col, 
                                   xlim, ylim,  
                                   main = "", sub = NULL, xlab,
                                   ylab, lty = par("lty"), lwd = par("lwd"), ...) {
  readpars <- par(no.readonly = TRUE)
  idx <- dimnames(x$lambda)
  lidx <- sapply(idx, length)
  names(lidx) <- names(idx) <- c("rows", "columns")

  dens <- array(NA, dim = lidx, dimnames = idx)
  for (i in idx[[2]]) {
    for (j in idx[[1]])
      dens[j,i] <- dnorm(x$lambda[j,i], mean = x$lambda[j,i],
                         sd = x$se[j,i])
  }

  if (missing(xlim)) {
    tmpL <- x$lambda - 3*x$se
    tmpH <- x$lambda + 3*x$se
    xlim <- c(min(tmpL), max(tmpH))
  } 
  
  if (by == "row") { 
    par(mfrow = c(lidx[1], 1), ...)
    
    if (missing(ylim)) { 
      ylim <- apply(dens, 1, max)
    }
    
    if (missing(xlab)) {
      xlab <- paste("Percentage", unlist(idx[[1]]))
    }
    else if (length(xlab) != length(idx[[1]])) {
      warning(paste("xlab needs to be", lidx[1], "long"))
      xlab <- rep(xlab, lidx[1])[1:lidx[1]]
    }
    if (missing(col)) {
      col <- rainbow(lidx[2])
    }
    else {
      if (length(col) < lidx[2])
        col <- rep(col, idx[[2]])[1:lidx[2]]
      if (length(col) > lidx[2])
        col <- col[1:lidx[2]]
    }
    names(col) <- idx[[2]]
    names(xlab) <- idx[[1]]
    
    for (ii in idx[[1]]) {
      x1 <- seq(tmpL[ii,1], tmpH[ii,1], by = 0.001)
      d1 <- dnorm(x1, mean = x$lambda[ii,1], sd = x$se[ii,1])
      plot(x1, d1, type = "l", xlim = xlim, col = col[1],
           ylim = c(0, max(ylim[ii])), main = main, sub = sub,
           xlab = xlab[ii], ylab = "Density", lty = lty,
           lwd = lwd)
      for (jj in idx[[2]][2:lidx[2]]) {
        x1 <- seq(tmpL[ii,jj], tmpH[ii,jj], by = 0.001)
        d1 <- dnorm(x1, mean = x$lambda[ii,jj], sd = x$se[ii,jj])
        lines(x1, d1, lty = lty, lwd = lwd, col = col[jj])
      }
      abline(v = 0, col = "grey50")
      abline(v = 1, col = "grey50")
    }
  }
  if (by == "column") {
    par(mfrow = c(lidx[2], 1), ...)

    if (missing(ylim)) { 
      ylim <- apply(dens, 2, max)
    }
    if (missing(xlab)) {
      xlab <- paste("Percentage", unlist(idx[[2]]))
    }
    else {
      if (length(xlab) != lidx[1]) {
        warning(paste("xlab needs to be", lidx[1], "long"))
        xlab <- rep(xlab, lidx[2])[1:lidx[2]]
      }
    }
    if (missing(col)) {
      col <- rainbow(lidx[1])
    }
    else {
      if (length(col) < lidx[1])
        col <- rep(col, idx[[1]])[1:lidx[1]]
      if (length(col) > lidx[1])
        col <- col[1:lidx[1]]
    }
    names(col) <- idx[[1]]
    names(xlab) <- idx[[2]]
    for (ii in idx[[2]]) {
      x1 <- seq(tmpL[1,ii], tmpH[1,ii], by = 0.001)
      d1 <- dnorm(x1, mean = x$lambda[1,ii], sd = x$se[1,ii])
      plot(x1, d1, type = "l", xlim = xlim, col = col[1],
           ylim = c(0, max(ylim[ii])), main = main, sub = sub,
           xlab = xlab[ii], ylab = "Density", lty = lty,
           lwd = lwd)
      for (jj in idx[[1]][2:lidx[1]]) {
        x1 <- seq(tmpL[jj,ii], tmpH[jj,ii], by = 0.001)
        d1 <- dnorm(x1, mean = x$lambda[jj,ii], sd = x$se[jj,ii])
        lines(x1, d1, lty = lty, lwd = lwd, col = col[jj])
      }
      abline(v = 0, col = "grey50")
      abline(v = 1, col = "grey50")
    }
  }
  return(invisible(x))
  par(readpars)
}
