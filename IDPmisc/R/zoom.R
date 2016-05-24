### zoom.R

zoom <- function(fun = plot, zoom.col="red", delay=3, ...){
  ##
  ## fun       plotting function
  ## zoom.col  color of clicked points
  ## delay     number of sec during which the 2nd zooming point is shown
  ##           before zoomed figure ist drawn
  ## ...       Arguments for plotting function
  ##
  ##
  ## Author: Rene Locher
  ## Version 03.06.05

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  options(locatorBell=FALSE)
  cat("Click mouse at corners of zoom area.\nRight click to stop zooming\n\n")

  fun(...)
  while(TRUE) {
    par(xpd=NA)
    p1 <- locator(n = 1, type="p", col=zoom.col)

    if(is.null(p1$x)){
      cat("\n")
      break
    }
    abline(v=p1$x,col=zoom.col)
    abline(h=p1$y,col=zoom.col)

    p2 <- locator(n = 1, type="p", col=zoom.col)
    if(is.null(p2$x)){
      cat("\n")
      break
    }
    abline(v=p2$x,col=zoom.col)
    abline(h=p2$y,col=zoom.col)
    Sys.sleep(delay)

    xx <- sort(c(p1$x,p2$x))
    yy <-  sort(c(p1$y,p2$y))

    ## check if zooming out
    usr <- par()$usr

    ## zooming out in x direction
    if (xx[1]<usr[1]) xx[1] <- usr[1]-diff(range(usr[1:2]/3))
    if (xx[2]<usr[1]) xx[2] <- usr[1]-diff(range(usr[1:2]/3))
    if (xx[1]>usr[2]) xx[1] <- usr[2]+diff(range(usr[1:2]/3))
    if (xx[2]>usr[2]) xx[2] <- usr[2]+diff(range(usr[1:2]/3))

    ## zooming out in y direction
    if (yy[1]<usr[3]) yy[1] <- usr[3]-diff(range(usr[3:4]/3))
    if (yy[2]<usr[3]) yy[2] <- usr[3]-diff(range(usr[3:4]/3))
    if (yy[1]>usr[4]) yy[1] <- usr[4]+diff(range(usr[3:4]/3))
    if (yy[2]>usr[4]) yy[2] <- usr[4]+diff(range(usr[3:4]/3))

    xlim <- range(xx)
    ylim <- range(yy)
    print(t(data.frame(xlim,ylim)))
    cat("\n")

    par(xpd=FALSE)
    fun(..., xlim = xlim, ylim = ylim)
  }

}  ## zoom

