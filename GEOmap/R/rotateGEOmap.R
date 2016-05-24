rotateGEOmap<-function(INmap, TARGlat,  TARGlon, LAT0,  LON0, beta=0)
  {
    ###   rotate a geomap list so the center is on the target lat lon
    ###     INmap   - geomap structure
    ###     TARGlat - target lat to rotate to 
    ###     TARGlon - target lon to rotate to 
    ###     LAT0    - center lat of map
    ###     LON0    - center lon of map

      ###   Beta (degrees) is a final rotation about the z axis

    if(missing(beta)) beta = 0

    
    v1 = ll2xyz(TARGlat,  TARGlon )
    v2 = ll2xyz(LAT0,  LON0)


    g = X.prod((v1), (v2))

    delta = (180/pi)*acos( sum(v1*v2)/(sqrt(sum(v1^2))*sqrt(sum(v2^2))))

    R1 =gmat(g, c(0,0,0) , -delta)

    if(beta != 0)
      {
       
        R2 =gmat(v1, c(0,0,0) , beta)

        
        R1 =  R1 %*% R2

      }


    MAP = INmap
    IN = 1:length(MAP$STROKES$index)
    for (i in IN) {
      j1 = MAP$STROKES$index[i] + 1
      j2 = j1 + MAP$STROKES$num[i] - 1
      if ((j1 > 0 & j2 > 0 & j2-j1 >= 0)) {
        JEC = j1:j2
      }
      xlon = MAP$POINTS$lon[JEC]
      ylat = MAP$POINTS$lat[JEC]
      
###   rotate these so that: KORLON0  KORLAT0 is the center
      ppp = Lll2xyz(ylat, xlon )
      
      X2 = cbind(ppp$x, ppp$y, ppp$z, rep(1, times=length(ppp$x)))
      
      XRot =      X2 %*% R1

      K = list(x= XRot[,1]  , y=XRot[,2],   z=XRot[,3])

      NEWLL = Lxyz2ll(K)

      MAP$POINTS$lon[JEC]=NEWLL$lon
      MAP$POINTS$lat[JEC]=NEWLL$lat
    }

    MAP = boundGEOmap(MAP)

    invisible(MAP)

  }

