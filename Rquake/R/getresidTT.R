
getresidTT<-function(Ldat,EQ, stas , vel)
{

   Ldat = LeftjustTime(Ldat)
     m1 = match(Ldat$name, stas$name)

        Ldat$lat = stas$lat[m1]
        Ldat$lon = stas$lon[m1]
        Ldat$z = stas$z[m1]

       
   MLAT = median(Ldat$lat)
   MLON = median(Ldat$lon)

   proj = GEOmap::setPROJ(type=2, LAT0=MLAT, LON0=MLON)

   XY = GEOmap::GLOB.XY(Ldat$lat, Ldat$lon, proj)
   Ldat$x = XY$x
   Ldat$y = XY$y

   eqxy = GEOmap::GLOB.XY(EQ$lat, EQ$lon, proj)
   EQ$x = eqxy$x
   EQ$y =  eqxy$y
   
   n = length(Ldat$sec)


    delx = EQ$x-Ldat$x
   dely = EQ$y-Ldat$y
     
   deltadis =sqrt( (delx)^2 +  (dely)^2)

   ## print(deltadis)

   RHS = rep(0, length = n)

   G1 = GETpsTT(Ldat$phase, eqz=EQ$z, staz=0, delx=delx, dely=dely,  deltadis=deltadis , vel)

     RHS = (Ldat$sec - EQ$sec - G1$TT)

   return(RHS)
 }
