##
## Plot errors in numerical inversion methods
##
## --------------------------------------------------------------------------

plot.rvgt.ierror <- function (x, maxonly=FALSE, tol=NA, ...)

  ## ------------------------------------------------------------------------
  ## Function for plotting graph of errors
  ## ------------------------------------------------------------------------
  ## x     : object of class "rvgt.ierror" containing quantiles of the
  ##         errors, or a list of such objects
  ## maxonly: if TRUE, only show maximal errors
  ## tol   : maximal tolerated error
  ## ------------------------------------------------------------------------
{
  ## check arguments
  if (!is.list(x))
    stop ("Invalid argument 'x'.")

  ## is this a number?
  if (!is.na(tol)) {
    tol <- as.numeric(tol)
    if (tol <= 0 )
      stop ("Invalid argument 'tol'.")
  }
  
  ## we have two cases:
  ## either result is an object of class "rvgt.ierror", or
  ## it is a list of such objects.
  ## The former case is transformed into the latter.
  if (class(x) == "rvgt.ierror") {
    ierr <- list(x)
  }
  else {
    ## here we should check the class of the list members
    if (class(x[[1]]) != "rvgt.ierror")
      stop ("Invalid argument 'x'.")
    ierr <- x
  }
  
  ## -- now create plot -----------------------------------------------------

  ## parameters for plot
  nplots <- length(ierr)

  ## y limits for ploting area
  err.max <- 0
  for (i in 1:nplots) {
    err.max <- max(err.max, ierr[[i]]$max)
  }

  if (!is.na(tol)) {
    err.max <- max(err.max, 1.1*tol)
    err.max <- min(err.max, 10*tol)
  }

  ## u limits for ploting area
  umin <- 1
  umax <- 0
  for (i in 1:nplots) {
    umin <- min(umin,ierr[[i]]$udomain[1])
    umax <- max(umax,ierr[[i]]$udomain[2])
  }

  ## create plotting aera with labels
  plot(1, err.max, xlim=c(umin,umax), ylim=c(0,err.max), type="n",
       xlab="u", ylab=ierr[[i]]$kind, ...)

  ## draw lines for each set of errors
  if (nplots>1) {
    ## more than one table --> draw maximal errors only 

    ## colors
    lc <- rainbow(nplots)

    for (i in 1:nplots) {
      res <- ierr[[i]]$res
      umin <- max(0,ierr[[i]]$udomain[1])
      umax <- min(1,ierr[[i]]$udomain[2])
      length <- (umax-umin) / res
      u <- umin + length * ((1:res)-0.5)

      lines(u, ierr[[i]]$max, type="l", col=lc[i])
      text(x=1, y=ierr[[i]]$med[res], col=lc[i],
           label=paste(" [",i,"]", sep=""))
    }
  }

  else {
    ## one table --> draw range for each interval

    ## resolution
    res <- ierr[[1]]$res

    ## domain
    umin <- max(0,ierr[[1]]$udomain[1])
    umax <- min(1,ierr[[1]]$udomain[2])

    ## u values
    length <- (umax-umin) / res
    u <- umin + length * ((1:res)-0.5)

    if (isTRUE(maxonly)) {
      lines(u, ierr[[1]]$max, type="l", col="red")
    }
    else {
      ## range
      polygon( c(u,rev(u)), c(ierr[[1]]$max, rev(ierr[[1]]$min)),
              col="#DFB5BD", border = NA )
      ## interquartile range
      polygon( c(u,rev(u)), c(ierr[[1]]$uqr, rev(ierr[[1]]$lqr)),
              col="#DC6E83", border = NA )
      ## median
      lines(u,ierr[[1]]$med,type="l",col="#CA1436")
      ## maximum
      lines(u,ierr[[1]]$max,type="l",col="red",lwd=1)
    }
  }

  ## line for given tolerance
  if (!is.na(tol))
    abline(tol,0,col="darkblue",lwd=2,lty=2)
}

## --------------------------------------------------------------------------

print.rvgt.ierror <- function (x, ...) {
  cat("\nrvgtest - ierror:\n")
  cat("   kind of error  =",x$kind,"\n")
  cat("   sample size    =",x$n,"\n");
  cat("   u-domain       = (",
      x$udomain[1],", ",x$udomain[2],")\n",sep="")
  cat("   max error      =",signif(max(x$max),4),"\n")
  cat("   mean abs error =",signif(sum(x$mad)/x$res,4),"\n")
  cat("   root of MSE    =",signif(sqrt(sum(x$mse)/x$res),4),"\n\n")
}

## --------------------------------------------------------------------------
