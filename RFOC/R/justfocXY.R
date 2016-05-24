`justfocXY` <-
function(MEC, x=x, y=y, size=c(1,1), fcol=gray(0.90), fcolback="white", xpd=TRUE)
  {
    if(missing(fcol)) { fcol = gray(0.90) }
    if(missing(fcolback)) { fcolback =  "white"  }
    if(missing(xpd)) {  xpd=TRUE}
    C = RPMG::circle()
    
    
    lines(x+size[1]*C$x,  y+size[2]*C$y, type='l', xpd=xpd)

    if(is.null(MEC$UP)) MEC$UP = TRUE
    ## pax = focpoint(MEC$P$az, MEC$P$dip,  lab="P", UP=MEC$UP, PLOT=FALSE)
    
    PLS = polyfoc(MEC$az1, MEC$dip1, MEC$az2, MEC$dip2, UP=MEC$UP, PLOT=FALSE)
    ## POK = list(x=PLS$Px, y =PLS$Py)
    ## kin = inpoly(pax$x, pax$y,POK)

    ## print(paste(sep=' ', kin, MEC$sense))
    
    if(MEC$sense==1)
      {
        polygon(x+size[1]*C$x,y+size[2]*C$y, col=fcolback, xpd=xpd )
        polygon(x+size[1]*PLS$Px, y+size[2]*PLS$Py, col=fcol, xpd=xpd )
      }else
    {
      polygon(x+size[1]*C$x,y+size[2]*C$y, col=fcol, xpd=xpd )
      polygon(x+size[1]*PLS$Px, y+size[2]*PLS$Py, col=fcolback, xpd=xpd )
    }

  }

