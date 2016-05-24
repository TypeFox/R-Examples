##
##  PURPOSE:   Plotting of computed univariate conditional cdf's
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   06/05/2010
##
##  FUNCTIONS: plot.NMixPredCondCDFMarg (06/05/2010)
##             
## ======================================================================

## *************************************************************
## plot.NMixPredCondCDFMarg
## *************************************************************
##
plot.NMixPredCondCDFMarg <- function(x, ixcond, imargin, prob, over=FALSE, auto.layout=TRUE, type="l", lwd=1, lty, col, qlwd=1, qlty, qcol, main, xlab, ylab, ylim, annot=TRUE, ...)
{  
  if (missing(ixcond) & missing(imargin)) stop("either ixcond or imargin must be given")

  input.miss.imargin <- missing(imargin)
  input.miss.ixcond  <- missing(ixcond)
  input.miss.ylim <- missing(ylim)
  
  p <- length(x$x)
  xcond <- x$x[[x$icond]]


  ##### Plot of one univariate conditional density given (in a sequel) selected (or all) values of the margin we have conditioned
  ##### -------------------------------------------------------------------------------------------------------------------------
  if (!missing(imargin)){
    if (length(imargin) != 1) stop("imargin must be a single number")
    if (imargin < 1 | imargin > p) stop(paste("imargin must be between 1 and ", p, sep=""))
    if (imargin == x$icond) stop(paste("imargin may not be equal to ", x$icond, sep=""))
   
    if (missing(ixcond)){
      ixcond <- 1:length(xcond)
    }else{  
      if (any(ixcond <= 0)) stop("ixcond must be all positive")
      if (any(ixcond > length(xcond))) stop(paste("too high value of ixcond, maximum is ", length(xcond), sep=""))
      xcond <- xcond[ixcond]
    }            

    nfig <- length(xcond)    
    if (auto.layout){
      if (over){
        oldPar <- par(mfrow=c(1, 1), bty="n")
      }else{
        LAY <- autolayout(nfig)
        oldPar <- par(bty="n")
        layout(LAY)
      }  
      on.exit(oldPar)
    }    

    if (missing(xlab)){
      xlab <- M1 <- paste("x", imargin, sep="")
      M2 <- paste("x", x$icond, sep="")
    }else{
      if (length(xlab) != 1 & length(xlab) != p) stop("incorrect xlab")
      if (length(xlab) == 1){
        M1 <- xlab
        M2 <- paste("x", x$icond, sep="")
      }  
      if (length(xlab) == p){
        M2 <- xlab[x$icond]
        xlab <- M1 <- xlab[imargin]
      }  
    }  
    #
    if (missing(ylab)) ylab <- "CDF"
    if (length(ylab) != 1 & length(ylab) != p) stop("incorrect ylab")
    if (length(ylab) == p) ylab <- ylab[imargin]
    #
    if (missing(main)){
      if (over) main <- paste(M1, " given ", M2, sep="")
      else      main <- paste(M1, " given ", M2, "=", round(xcond, 2), sep="")
    }  
    if (length(main) != 1 & length(main) != length(xcond)) stop("incorrect main")
    if (over){
      if (length(main) > 1) main <- main[1]
    }else{
      if (length(main) == 1) main <- rep(main, length(xcond))
    }  
    #
    if (missing(col)){
      if (over){
        col <- heat_hcl(length(xcond), h=c(0, -100), c.=c(40, 80), l=c(75, 40), power=1)             ### sequential palette
        col <- col[length(col):1]
      }  
      else      col <- rep("darkblue", length(xcond))
    }
    #
    if (missing(qcol)) qcol <- col
    #
    if (length(col) != 1 & length(col) != length(xcond)) stop("incorrect col")
    if (length(col) == 1) col <- rep(col, length(xcond))    
    #
    if (missing(lty)) lty <- rep(1, length(xcond))
    if (length(lty) != 1 & length(lty) != length(xcond)) stop("incorrect lty")
    if (length(lty) == 1) lty <- rep(lty, length(xcond))    
    #    
    if (length(lwd) != 1 & length(lwd) != length(xcond)) stop("incorrect lwd")
    if (length(lwd) == 1) lwd <- rep(lwd, length(xcond))    
    #
    if (length(qcol) != 1 & length(qcol) != length(xcond)) stop("incorrect qcol")
    if (length(qcol) == 1) qcol <- rep(qcol, length(xcond))    
    #
    if (missing(qlty)) qlty <- rep(2, length(xcond))
    if (length(qlty) != 1 & length(qlty) != length(xcond)) stop("incorrect qlty")
    if (length(qlty) == 1) qlty <- rep(qlty, length(xcond))    
    #    
    if (length(qlwd) != 1 & length(qlwd) != length(xcond)) stop("incorrect qlwd")
    if (length(qlwd) == 1) qlwd <- rep(qlwd, length(xcond))        
    
    if (over){
      yy <- x$cdf[[ixcond[1]]][[imargin]]
      if (length(xcond) > 1) for (i in 2:length(xcond)) yy <- c(yy, x$cdf[[ixcond[i]]][[imargin]])
      if (input.miss.ylim) ylim <- c(0, 1)
      plot(x$x[[imargin]], x$cdf[[ixcond[1]]][[imargin]], type=type, col=col[1], lty=lty[1], lwd=lwd[1], main=main, xlab=xlab, ylab=ylab, ylim=ylim)
      if (!missing(prob)){
        for (j in 1:length(prob)){
          qnaam <- paste("q", prob[j]*100, "%", sep="")
          if (qnaam %in% names(x)) lines(x$x[[imargin]], x[[qnaam]][[ixcond[1]]][[imargin]], col=qcol[1], lty=qlty[1], lwd=qlwd[1])
        }  
      }
      
      if (length(xcond) > 1) for (i in 2:length(xcond)){
        lines(x$x[[imargin]], x$cdf[[ixcond[i]]][[imargin]], col=col[i], lty=lty[i], lwd=lwd[i])
        if (!missing(prob)){
          for (j in 1:length(prob)){
            qnaam <- paste("q", prob[j]*100, "%", sep="")
            if (qnaam %in% names(x)) lines(x$x[[imargin]], x[[qnaam]][[ixcond[i]]][[imargin]], col=qcol[i], lty=qlty[i], lwd=qlwd[i])
          }  
        }        
      }
      if (annot){
        LEG <- paste(M2, "=", round(xcond, 2))
        legend(x$x[[imargin]][1], ylim[2], legend=LEG, lty=lty, col=col, bty="n")
      }  
    }else{
      for (i in 1:length(xcond)){
        if (input.miss.ylim) ylim <- c(0, 1)
        plot(x$x[[imargin]], x$cdf[[ixcond[i]]][[imargin]], type=type, col=col[i], lty=lty[i], lwd=lwd[i], main=main[i], xlab=xlab, ylab=ylab, ylim=ylim)
        if (!missing(prob)){
          for (j in 1:length(prob)){
            qnaam <- paste("q", prob[j]*100, "%", sep="")
            if (qnaam %in% names(x)) lines(x$x[[imargin]], x[[qnaam]][[ixcond[i]]][[imargin]], col=qcol[i], lty=qlty[i], lwd=qlwd[i])
          }  
        }        
      }  
    }  
  }  
    
  
  ##### Plot of all univariate conditional densities given one value of the margin we have conditioned
  ##### ------------------------------------------------------------------------------------------------
  if (input.miss.imargin & !input.miss.ixcond){
    if (length(ixcond) != 1) stop("ixcond must be a single number")
    if (ixcond < 1 | ixcond > length(xcond)) stop(paste("ixcond must be between 1 and ", length(xcond), sep=""))
    
    nfig <- p - 1
    if (auto.layout){    
      LAY <- autolayout(nfig)    
      oldPar <- par(bty="n")
      layout(LAY)
      on.exit(oldPar)
    }
    
    xcond <- xcond[ixcond]       ### we will plot densities X[d] | X[icond] = xcond
    #
    if (missing(xlab)) xlab <- paste("x", 1:p, sep="")
    if (length(xlab) != 1 & length(xlab) != p) stop("incorrect xlab")    
    if (length(xlab) == 1) xlab <- rep(xlab, p)
    #
    if (missing(ylab)) ylab <- "CDF"
    if (length(ylab) != 1 & length(ylab) != p) stop("incorrect ylab")    
    if (length(ylab) == 1) ylab <- rep(ylab, p)
    #
    if (missing(main)) main <- paste(xlab, " given ", xlab[x$icond], "=", round(xcond, 2), sep="")
    if (length(main) != 1 & length(main) != p) stop("incorrect main")        
    if (length(main) == 1) main <- rep(main, p)
    #
    if (missing(col)) col <- "darkblue"
    if (length(col) != 1 & length(col) != p) stop("incorrect col")
    if (length(col) == 1) col <- rep(col, p)    
    #
    if (missing(lty)) lty <- 1
    if (length(lty) != 1 & length(lty) != p) stop("incorrect lty")
    if (length(lty) == 1) lty <- rep(lty, p)
    #
    if (length(lwd) != 1 & length(lwd) != p) stop("incorrect lwd")    
    if (length(lwd) == 1) lwd <- rep(lwd, p)    
    #
    if (missing(qcol)) qcol <- col
    if (length(qcol) != 1 & length(qcol) != p) stop("incorrect qcol")
    if (length(qcol) == 1) qcol <- rep(qcol, p)    
    #
    if (missing(qlty)) qlty <- 1
    if (length(qlty) != 1 & length(qlty) != p) stop("incorrect qlty")
    if (length(qlty) == 1) qlty <- rep(qlty, p)
    #
    if (length(qlwd) != 1 & length(qlwd) != p) stop("incorrect qlwd")    
    if (length(qlwd) == 1) qlwd <- rep(qlwd, p)        

    for (i in (1:p)[-x$icond]){
      plot(x$x[[i]], x$cdf[[ixcond]][[i]], type=type, col=col[i], lty=lty[i], lwd=lwd[i], main=main[i], xlab=xlab[i], ylab=ylab[i])
      if (!missing(prob)){
        for (j in 1:length(prob)){
          qnaam <- paste("q", prob[j]*100, "%", sep="")
          if (qnaam %in% names(x)) lines(x$x[[i]], x[[qnaam]][[ixcond]][[i]], col=qcol[i], lty=qlty[i], lwd=qlwd[i])
        }  
      }      
    }      
  }  

  return(invisible(x))    
}

