`plotUTM` <-
function(proj, LIM, shiftlon=0)
{

  if(missing(shiftlon)) { shiftlon=0 }

  utmx = seq(from=(-180), to=180, by=6)
  utmy = seq(from =-56, to =72, by=8)

  LONS = RPMG::fmod(utmx, 360)

  htags =  LETTERS[6:23]
  htags = htags[-c(4, 10)]


  lon180 = RPMG::fmod(c(LIM[1], LIM[3]), 360 )

lon180[lon180>180] = lon180[lon180>180]-360
  
  
 ##### ulats = utmy[utmy>=LIM[2] & utmy<=LIM[4] ]
 ##### ulons = utmx[LONS>=LIM[1] & LONS<=LIM[3] ]

  ux1 = findInterval(lon180[1], utmx)
  ux2 = findInterval(lon180[2], utmx)
   ulons =utmx[ux1:ux2]

  uy1 = findInterval(LIM[2], utmy)
  uy2 = findInterval(LIM[4], utmy)
  
  ulats =utmy[uy1:uy2]


  j = RPMG::fmod(seq(from=(-180), to=180, by=6)-shiftlon, 360)
  jhalf = RPMG::fmod(seq(from=(-177), to=177, by=6)-shiftlon, 360)
###  text(jhalf , u[3], labels=1:length(jhalf), pos=3)


  LONLABS =1:length(jhalf) 
  
  for(i in 1:length(ulons))
    {
      XY = GLOB.XY(ulats, rep(ulons[i], length(ulats)) , proj)
      lines(XY, col='blue')
      

    }
  
  for(i in 1:length(ulats))
    {
      XY = GLOB.XY(rep(ulats[i], length( ulons))  , ulons , proj)
      lines(XY, col='blue')

    }


  halflat = utmy[1:(length(utmy)-1)]+diff(utmy)/2
  
  XY = GLOB.XY( halflat[uy1:uy2] ,  rep(ulons[1], length(ulats)), proj)
text(XY, labels=htags[uy1:uy2])

  

  XY = GLOB.XY(  rep(ulats[i], length( ulons)) , jhalf[ux1:ux2], proj)
text(XY, labels=LONLABS[ux1:ux2])
  
######  pusa = plotusa(); proj = pusa$PROJ; LIM=pusa$LIM
######   plotUTM(proj, LIM)
######   proj = pusa$PROJ

}

