`GETXprofile` <-
function(jx, jy, jz, LAB="A", myloc=NULL, PLOT=FALSE, NEWDEV=TRUE,  asp=1)
{
###############  get a cross section through a topographic DEM
############   GETXprofile : xsec  through topo

  if(missing(myloc)) {  myloc=NULL }
    if(missing(LAB)) {   LAB="A" }
    if(missing(asp)) {   asp=NULL }

  if(is.null(myloc))
    {
      myloc = locator(2, type='o')
    }

  delx = mean(abs(diff(jx)))
  dely = mean(abs(diff(jy)))

MDX = diff(myloc$x)
MDY = diff(myloc$y)
  
  dx = sign(MDX)*delx
  dy = sign(MDY)*dely

  
  ####  this is overkill - need simply slope-intercept
  
  ###Lrunvent = lm ( myloc$y ~ myloc$x )

if(myloc$x[1]==myloc$x[2])
  {


    newy =  seq(from=myloc$y[1], to=myloc$y[2], by=dy)
    JX = rep(myloc$x[1], length=length(newy))
      allx = JX
  ally = newy

  ox = order(ally)
  allx = allx[ox]
  ally = ally[ox]
  
  }
  else
    {
  slope = MDY/MDX
  intercept = myloc$y[1] - slope*myloc$x[1]

  newx = seq(from=myloc$x[1], to=myloc$x[2], by=dx)
  
  JY = intercept+slope*newx

  newy =  seq(from=myloc$y[1], to=myloc$y[2], by=dy)

  JX = (newy-intercept)/slope


  allx = c(newx, JX)
  ally = c(JY, newy)

  ox = order(allx)
  allx = allx[ox]
  ally = ally[ox]
  }
###  points(jx, JY)
#####  here we could use findIntervl -
  ###  but it is overkill since the intervals are
  ###  uniform.  So - go with simpler approach
  
  iboxx = findInterval(allx, jx)

  iboxy =  findInterval(ally, jy, all.inside = FALSE)


#####



  flag = (iboxy>0 & iboxy<length(jy))  &  (iboxx>0 & iboxx<length(jx)) 
 
  

### points(jx[flag],JY[flag], col=p2)

  pts = cbind(iboxx[flag], iboxy[flag])

  LX = allx[flag]
  LY = ally[flag]

  RX = sqrt((LX-LX[1])^2+ (LY-LY[1])^2)

  LZ = jz[pts]

  pnt1 = sqrt((myloc$x[1]-LX[1])^2+    (myloc$y[1]-LY[1])^2)
  pnt2 = sqrt((myloc$x[2]-LX[1])^2+    (myloc$y[2]-LY[1])^2)


  ###  I do not recall why this is in here 
 ## px1 = findInterval(pnt1, RX)
 ## px2 = findInterval(pnt2, RX)

  px1 = which.min(abs(RX-pnt1))
  px2 =  which.min(abs(RX-pnt2))


  dist = sqrt((myloc$x[2]-myloc$x[1])^2+(myloc$y[2]-myloc$y[1])^2)

BX =  myloc$x[1]  +  RX*(myloc$x[2]-myloc$x[1])/dist


BY =  myloc$y[1]  + RX*(myloc$y[2]-myloc$y[1])/dist



  if(PLOT==TRUE)
    {
      ###screens(2)
      if(NEWDEV) cdev = dev.cur()
       if(NEWDEV)  dev.new()
      plot(RX, LZ, type='l', xlab="m", ylab="m", ylim=range(jz, na.rm=TRUE), asp=asp)
      points( c(pnt1, pnt2), c(LZ[px1], LZ[px2] ), pch=c(6,8), col=c(2,4) )

      text(c(pnt1, pnt2),c(LZ[px1], LZ[px2] ) , labels=c(LAB, paste(sep="", LAB, "'")), col=c(2,4), pos=3 )
      if(NEWDEV)  dev.set(cdev)
    }
  invisible(list(RX=RX, RZ=LZ, LOC=myloc, LOCx=BX, LOCy=BY, LAB=LAB ))

}

