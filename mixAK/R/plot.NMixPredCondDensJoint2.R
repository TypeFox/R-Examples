##
##  PURPOSE:   Plotting of computed predictive pairwise bivariate conditional densities
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   01/06/2009
##  LOG:       08/11/2011:  add.contour, col.add.contour arguments added
##
##  FUNCTIONS: plot.NMixPredCondDensJoint2 (01/06/2009)
##             
## ======================================================================

## *************************************************************
## plot.NMixPredCondDensJoint2
## *************************************************************
##
plot.NMixPredCondDensJoint2 <- function(x, ixcond, imargin, contour=FALSE, add.contour=TRUE, col.add.contour="brown", auto.layout=TRUE, col, lwd=1, main, xylab, ...)
{
  if (missing(ixcond) & missing(imargin)) stop("either ixcond or imargin must be given")

  input.miss.imargin <- missing(imargin)
  input.miss.ixcond  <- missing(ixcond)
  
  p <- length(x$x)  
  xcond <- x$x[[x$icond]]

  miss.main <- missing(main)
  miss.xylab <- missing(xylab)
  if (!miss.xylab){
    if (length(xylab) != p) stop("xylab must be of length", p)
  }  
  
  if (missing(col)){
    if (contour) col <- "darkblue"
    else{
      col <- rev(heat_hcl(33, c.=c(80, 30), l=c(30, 90), power=c(1/5, 1.3)))
    }  
  }


  ##### Plot of one bivariate conditional density given (in a sequel) selected (or all) values of the margin we have conditioned
  ##### ------------------------------------------------------------------------------------------------------------------------
  if (!missing(imargin)){
    if (length(imargin) != 2) stop("imargin must be a vector of length 2")
    if (any(imargin < 1) | any(imargin > p)) stop(paste("imargin must be both between 1 and ", p, sep=""))
    if (any(imargin == x$icond)) stop(paste("imargin may not be equal to ", x$icond, sep=""))
    if (imargin[1] == imargin[2]) stop("imargin must be two different numbers")
    
    if (imargin[1] > imargin[2]) NAAM <- paste(imargin[2], "-", imargin[1], sep="")
    else                         NAAM <- paste(imargin[1], "-", imargin[2], sep="")

    if (miss.xylab){
      xlab <- paste("x", imargin[1], sep="")
      ylab <- paste("x", imargin[2], sep="")
      zlab <- paste("x", x$icond, sep="") 
    }else{
      xlab <- xylab[imargin[1]]
      ylab <- xylab[imargin[2]]
      zlab <- xylab[x$icond]
    }

    if (missing(ixcond)){
      ixcond <- 1:length(xcond)
    }else{  
      if (any(ixcond <= 0)) stop("ixcond must be all positive")
      if (any(ixcond > length(xcond))) stop(paste("too high value of ixcond, maximum is ", length(xcond), sep=""))
      xcond <- xcond[ixcond]
    }            
    
    nfig <- length(xcond)
    if (auto.layout){
      LAY <- autolayout(nfig)
      oldPar <- par(bty="n")
      layout(LAY)
      on.exit(oldPar)
    }    

    if (miss.main){
      main <- paste("(", xlab, ", ", ylab, ") given ", zlab, "=", round(xcond, 2), sep="")
    }
    if (length(main) == 1) main <- rep(main, length(xcond))
    if (length(main) != length(xcond)) stop(paste("main must be of length ", length(xcond), sep=""))
    
    for (i in 1:length(xcond)){
      if (imargin[1] > imargin[2]) dx <- t(x$dens[[ixcond[i]]][[NAAM]])
      else                         dx <- x$dens[[ixcond[i]]][[NAAM]]

      if (contour){
        contour(x$x[[imargin[1]]], x$x[[imargin[2]]], dx, col=col, main=main[i], xlab=xlab, ylab=ylab, lwd=lwd, ...)
      }else{
        image(x$x[[imargin[1]]], x$x[[imargin[2]]], dx, col=col, main=main[i], xlab=xlab, ylab=ylab, ...)
        if (add.contour) contour(x$x[[imargin[1]]], x$x[[imargin[2]]], dx, col=col.add.contour, lwd=lwd, add=TRUE)
      }                        
    }  
  }

  
  ##### Plot of all bivariate conditional densities given one value of the margin we have conditioned
  ##### ------------------------------------------------------------------------------------------------
  if (input.miss.imargin & !input.miss.ixcond){  
    if (length(ixcond) != 1) stop("ixcond must be a single number")
    if (ixcond < 1 | ixcond > length(xcond)) stop(paste("ixcond must be between 1 and ", length(xcond), sep=""))
    
    nfig <- (p - 1)*(p - 2)/2
    if (auto.layout){    
      LAY <- autolayout(nfig)    
      oldPar <- par(bty="n")
      layout(LAY)
      on.exit(oldPar)
    }

    xcond <- xcond[ixcond]       ### we will plot densities X[d] | X[icond] = xcond

    for (m0 in 1:(p-1)){
      for (m1 in (m0+1):p){
        if (m0 == x$icond | m1 == x$icond) next

        NAAM <- paste(m0, "-", m1, sep="")  
        if (miss.xylab){
          xlab <- paste("x", m0, sep="")
          ylab <- paste("x", m1, sep="")
          zlab <- paste("x", x$icond, sep="") 
        }else{
          xlab <- xylab[m0]
          ylab <- xylab[m1]
          zlab <- xylab[x$icond]
        }
        if (miss.main){
          main <- paste("(", xlab, ", ", ylab, ") given ", zlab, "=", round(xcond, 2), sep="")
        }
        
        if (contour){
          contour(x$x[[m0]], x$x[[m1]], x$dens[[ixcond]][[NAAM]], col=col, main=main, xlab=xlab, ylab=ylab, lwd=lwd, ...)
        }else{
          image(x$x[[m0]], x$x[[m1]], x$dens[[ixcond]][[NAAM]], col=col, main=main, xlab=xlab, ylab=ylab, ...)
          if (add.contour) contour(x$x[[m0]], x$x[[m1]], x$dens[[ixcond]][[NAAM]], col=col.add.contour, lwd=lwd, add=TRUE)
        }                  
      }
    }          
  }  
     
  return(invisible(x)) 
}  
