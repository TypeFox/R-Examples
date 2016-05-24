explode<-function(fxy, dixplo=1, mult=1, cenx=0, ceny=0, PLOT=FALSE)
  {
############   calculate the location of an exploded set of points
    ### cenx and ceny as the center of the explosion,
    ###  if missing use the mean of the coordinates
    
    if(missing(mult)) { mult = 1 }
    if(missing(cenx)) {cenx=mean(fxy$x) }
   if(missing(ceny)) {ceny=mean(fxy$y) }
   if(missing(PLOT)) {PLOT=FALSE }

    

    
    mfxy = list(x=cenx,y=ceny)

    dis1 = sqrt( (fxy$x-mfxy$x)^2 + (fxy$y-mfxy$y)^2)

    if(missing(dixplo) ) {  dixplo  = mean(dis1) }
    dixplo  = mult * dixplo 
    DX = (fxy$x-mfxy$x)/dis1
    DY = (fxy$y-mfxy$y)/dis1

    PX = fxy$x+ dixplo* DX
    PY = fxy$y+ dixplo* DY

    if(PLOT)
      {
        plot(range(fxy$x, PX) , range(PY,  fxy$y), asp=1, type='n')
        points(fxy$x, fxy$y , col='green', pch=3)
        points(PX, PY, col='red', pch=6,  xpd=TRUE)

        segments(fxy$x, fxy$y,PX, PY)
        points(cenx, ceny, pch=8, col='black')

        
      }

    return(list(x=PX, y=PY))
  }

