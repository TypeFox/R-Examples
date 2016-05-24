##==============================================================================
## colorlegend  : adds a color legend to a plot
##==============================================================================

colorlegend <- function(col=femmecol(100), zlim, zlevels=5,
  dz=NULL, zval=NULL, log=FALSE, posx=c(0.9,0.93), posy=c(0.05,0.9),
  main=NULL, main.cex=1.0, main.col="black", lab.col="black",
  digit=0, left=FALSE, ...) {

  ncol   <- length (col)
  par (new=TRUE)
  omar <- nmar <- par("mar")
  nmar[c(2,4)]<-0
  par (mar = nmar)

  emptyplot()
  pars   <- par("usr")

  ## Rectangle positions on x and y-axis
  dx     <- pars[2]-pars[1]
  xmin   <- pars[1]+posx[1]*dx
  xmax   <- pars[1]+posx[2]*dx

  dy     <- pars[4]-pars[3]
  ymin   <- pars[3]+posy[1]*dy
  ymax   <- pars[3]+posy[2]*dy

  ## z-values
  if (!is.null(zval)) {
    zz<-zval
    dz<-NULL
  }

  if (is.null(dz)&is.null(zval))
    if (! is.null(zlevels)) {
      if (log) {
        zz <- 10^(pretty(log10(zlim),n=(zlevels+1)))
      } else
        zz <-     pretty(zlim,n=(zlevels+1))
    } else zz <- NULL
  if (!is.null(dz)) {
    if (log)
      zz <- 10^(seq(log10(zlim[1]),log10(zlim[2]),by=dz))
    if (!log)
      zz <- seq(zlim[1],zlim[2],by=dz)
  }

  if (log) {
    zlim <- log10(zlim)
    if (! is.null(zz))
      zz   <- log10(zz)
  }

  zmin   <- zlim[1]
  zmax   <- zlim[2]

  ## colors
  Y <- seq(ymin,ymax,length.out=ncol+1)
  rect(xmin,Y[-(ncol+1)],xmax,Y[-1],col=col,border=NA)
  rect(xmin,ymin,xmax,ymax,border=lab.col)

  if (! is.null(zz)) {
  ## labels
    dx     <- (xmax-xmin)
    dy     <- (ymax-ymin)

    if (left) {
      Dx  <-  -dx  # labels on left..
      pos <-   2
      xpos <- xmin+Dx*0.5
    } else {
      Dx  <- +dx  # labels on right..
      pos <- 4
      xpos <- xmax+Dx*0.5
    }

    Ypos <- ymin+(zz-zmin)/(zmax-zmin)*dy
    segments(xmin,Ypos,xmax,Ypos,col=lab.col)
    segments(xpos+Dx*0.25,Ypos,xmin,Ypos,col=lab.col)
    text (xpos,Ypos,formatC(zz,digits=digit,format="f"),pos=pos,col=lab.col,...)
  }

  if  (!is.null(main)) {
    for (i in length(main):1)
       text (x=mean(c(xmin,xmax)),y=ymax+0.05*(length(main)-i+1),
             labels=main[i],adj=c(0.5,0.5),cex=main.cex,col=main.col)
  }
  par (new=FALSE)
  par (mar=omar)

}
