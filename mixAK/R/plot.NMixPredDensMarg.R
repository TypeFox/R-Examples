##
##  PURPOSE:   Plotting of computed predictive marginal (univariate) densities
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   03/12/2007
##             02/04/2015:  covariates on mixture weights allowed
##
##  FUNCTIONS: plot.NMixPredDensMarg (03/12/2007)
##             
## ======================================================================

## *************************************************************
## plot.NMixPredDensMarg
## *************************************************************
plot.NMixPredDensMarg <- function(x, K=0, auto.layout=TRUE, type="l", col="darkblue", lty=1, lwd=1, main, xlab, ylab, ...)
{
  if (is.null(x$nx_w)) x$nx_w <- 1
  if (is.null(x$lx_w)) x$lx_w <- ""
  
  p <- length(x$x) / x$nx_w
  nfig <- length(x$x)

  miss.main <- missing(main)
  
  if (missing(xlab)){
    if (p == 1) xlab <- "x"
    else        xlab <- paste("x", 1:p, sep="")
  }else{
    if (length(xlab) == 1){
      xlab <- rep(xlab, p)
    }else{
      if (length(xlab) != p) stop(paste("xlab must be of length 1 or ", p, sep=""))
    }  
  }  
  
  if (missing(ylab)){
    ylab <- rep("Density", p)
  }else{
    if (length(ylab) == 1){
      ylab <- rep(ylab, p)
    }else{
      if (length(ylab) != p) stop(paste("ylab must be of length 1 or ", p, sep=""))
    }  
  }  

  #if (length(K) != 1) stop("K must be of length 1")
  if (any(K < 0)) stop("K must not be negative")
  if (any(K > length(x$densK[[1]]))) stop("K is too high")

  percK <- paste(" (", round(x$propK*100, 1), "%)", sep="")
      
  if (auto.layout){    
    LAY <- autolayout(nfig)    
    oldPar <- par(bty="n")
    layout(LAY)
    on.exit(oldPar)
  }

  if (length(K) == 1){
    if (x$nx_w == 1){
      for (m0 in 1:p){
        if (K == 0){
          dx <- x$dens[[paste(m0)]]
          #main2  <- if (is.null(x$MCMC.length)) "" else paste(" (MCMC length = ", x$MCMC.length, ")", sep="")
          #main2b <- if (is.null(x$MCMC.length)) "" else paste("MCMC length = ", x$MCMC.length, sep="")
          main2 <- ""
          main2b <- ""
        }else{
          dx <- x$densK[[paste(m0)]][[K]]
          main2  <- paste(",  K = ", K, percK[K], sep="")
          main2b <- paste("K = ", K, percK[K], sep="")                  
        }
      
        if (miss.main){
          if (p == 1) MAIN <- main2b
          else        MAIN <- paste("Margin ", m0, main2, sep="")
        }else{
          if (length(main) == 1) MAIN <- main
          else{
            if (length(main) != p) stop(paste("main must be of length 1 or ", p, sep=""))
            MAIN <- main[m0]
          }  
        }  

        plot(x$x[[m0]], dx, type=type, col=col[1], main=MAIN, xlab=xlab[m0], ylab=ylab[m0], lty=lty[1], lwd=lwd[1], ...)
      }
    }else{     ## x$nx_w > 1
      for (ixw in 1:x$nx_w){
        for (m0 in 1:p){
          if (K == 0){
            dx <- x$dens[[paste(m0, "-", x$lx_w[ixw], sep = "")]]
            #main2  <- if (is.null(x$MCMC.length)) "" else paste(" (MCMC length = ", x$MCMC.length, ")", sep="")
            #main2b <- if (is.null(x$MCMC.length)) "" else paste("MCMC length = ", x$MCMC.length, sep="")
            main2 <- ""
            main2b <- ""
          }else{
            dx <- x$densK[[paste(m0, "-", x$lx_w[ixw], sep = "")]][[K]]
            main2  <- paste(",  K = ", K, percK[K], sep="")
            main2b <- paste("K = ", K, percK[K], sep="")                  
          }
      
          if (miss.main){
            if (p == 1) MAIN <- paste(x$lx_w[ixw], main2b, sep = "")
            else        MAIN <- paste(x$lx_w[ixw], ": Margin ", m0, main2, sep="")
          }else{
            if (length(main) == 1) MAIN <- main
            else{
              if (length(main) != p) stop(paste("main must be of length 1 or ", p, sep=""))
              MAIN <- main[m0]
            }  
          }  

          plot(x$x[[paste("x", m0, "-", x$lx_w[ixw], sep = "")]], dx, type=type, col=col[1], main=MAIN, xlab=xlab[m0], ylab=ylab[m0], lty=lty[1], lwd=lwd[1], ...)
        }
      }  
    }          ## end of x$nx_w > 1

  }else{     ### length(K) > 1
    if (length(col) == 1) col <- rep(col, length(K))
    if (length(col) != length(K)) stop("incorrect length(col) argument")
    if (length(lty) == 1) lty <- rep(lty, length(K))
    if (length(lty) != length(K)) stop("incorrect length(lty) argument")
    if (length(lwd) == 1) lwd <- rep(lwd, length(K))
    if (length(lwd) != length(K)) stop("incorrect length(lwd) argument")

    if (x$nx_w == 1){
      for (m0 in 1:p){
        if (miss.main){
          if (p == 1) main <- if (is.null(x$MCMC.length)) "" else paste("MCMC length = ", x$MCMC.length, sep="")
          else        main <- if (is.null(x$MCMC.length)) "" else paste("Margin ", m0, " (MCMC length = ", x$MCMC.length, ")", sep="")
        }else{
          if (length(main) == 1) MAIN <- main
          else{
            if (length(main) != p) stop(paste("main must be of length 1 or ", p, sep=""))
            MAIN <- main[m0]
          }
        }  
      
        for (kk in 1:length(K)){
          if (K[kk] == 0) dx <- x$dens[[paste(m0)]]
          else            dx <- x$densK[[paste(m0)]][[K[kk]]]
          if (kk == 1) plot(x$x[[m0]], dx, type=type, col=col[1], main=main, xlab=xlab[m0], ylab=ylab[m0], lty=lty[1], lwd=lwd[1], ...)
          else         lines(x$x[[m0]], dx, col=col[kk], lty=lty[kk], lwd=lwd[kk])
        } 
      }
    }else{     ## x$nx_w > 1
      stop("not implemented for varying number of K's")
    }          ## end of x$nx_w > 1
  }  
    
  return(invisible(x))   
}  


