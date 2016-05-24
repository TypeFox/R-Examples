EQXYresid <-
function(XY, vel=list() , h1=c(0,0,0,0)  , PLOT=FALSE)
  {
###  init is a list with lat lon z and sec (sec is relative to the minute reference)
    if(missing(PLOT))  { PLOT=TRUE }
   
    if(missing(vel))
      {
       
    LITHOS.vel=defaultVEL(2)
       
        vel= LITHOS.vel
      }


    N = length(XY$phase)


    delx = XY$x-h1[1]
    dely = XY$y-h1[2]
    XY$r = sqrt((delx)^2 + (dely)^2)


    
  ###  XY$c = rep(0, length(XY$r))
  ###  XY$s = rep(0, length(XY$r))

   ### XY$c[XY$r>0] = (XY$x[XY$r>0]-h1[1])/XY$r[XY$r>0]
   ### XY$s[XY$r>0] = (XY$y[XY$r>0]-h1[2])/XY$r[XY$r>0]

  G1 = GETpsTT(XY$phase, eqz=h1[3], staz=0, delx=delx, dely=dely,  deltadis=XY$r , vel)

    rhs = rep(NA, times=N)

    rhs = XY$sec- (h1[4] + G1$TT+XY$cor ) 


  return( rhs)

  }
