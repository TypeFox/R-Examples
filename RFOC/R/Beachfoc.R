`Beachfoc` <-
function(MEC, fcol=gray(0.90), fcolback="white", ALIM = c(-1, -1, +1, +1))
{

  if(missing(fcol)) { fcol = gray(0.90) }
  if(missing(fcolback)) { fcolback =  "white"  }
  if(missing(ALIM))  {ALIM=c(-1,-1, +1, +1) }

  
  if(is.null(MEC$UP)) MEC$UP=TRUE
  if(is.null(MEC$LIM)) MEC$LIM = ALIM
  
  
##### net(1)
  C = RPMG::circle()
 #####  print(MEC$LIM)
  
  plot(C$x,C$y, type='n', axes=FALSE, asp=1, xlab="", ylab="", xlim=c(MEC$LIM[1],MEC$LIM[3]), ylim=c(MEC$LIM[2],MEC$LIM[4]) )
  lines(C$x,C$y, type='l')
  
####circ.tics()

 
  pax = focpoint(MEC$P$az, MEC$P$dip,  lab="P", UP=MEC$UP, PLOT=FALSE)
  
  PLS = polyfoc(MEC$az1, MEC$dip1, MEC$az2, MEC$dip2, UP=MEC$UP)
  
###   POK = list(x=PLS$Px, y =PLS$Py)
  
  ###  OLD, does not work: kin = inpoly(pax$x, pax$y,POK)
 kin = inout(cbind(pax$x, pax$y) ,cbind(PLS$Px, y =PLS$Py), bound=TRUE)

  if(kin==0)
    {
      polygon(PLS$Px, PLS$Py, col=fcol )
    }
  else
  {
    polygon(C$x,C$y, col=fcol )
    polygon(PLS$Px, PLS$Py, col=fcolback )
  }

###  lines(PLS$Px, PLS$Py, col="black", lwd=2)

  
}

