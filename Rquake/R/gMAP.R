gMAP <-
function(nh, g, idev=3)
  {

    
    
    Apf =  nh$pickfile

    if(is.null(Apf))
      {
        print("gMAP: no pickfile")
        dev.set( g$MAINdev)
        return()
      }

    sta = nh$sta

    MASTA  = match(Apf$STAS$name, sta$name)
    
    stalats =  sta$lat[MASTA]
    stalons = sta$lon[MASTA]
    stanam  = sta$name[MASTA]

  ###  srclat = Apf$LOC$lat
  ###  srclon = Apf$LOC$lon

    astas = nh$STNS
    acols = g$pcols[nh$pcol]

    colmatch = match(stanam, astas)
    
    pcols  = g$pcols[colmatch]
    
    
    print("DOING MAP")


    dev.set( dev.next() )

    print(sta)
    
    gproj = doAmap(sta)
    
    contPFarrivals(Apf, sta, proj=gproj, image=TRUE , add=TRUE)
  if(!is.na(Apf$LOC$lat))
    {
gstaxy = GEOmap::GLOB.XY(stalats, stalons, gproj)
gevxy = GEOmap::GLOB.XY(Apf$LOC$lat, Apf$LOC$lon, gproj)
 points(gevxy$x, gevxy$y,  pch=8, col='blue')
 segments(gevxy$x, gevxy$y, gstaxy$x, gstaxy$y, col='red')
    }
  
 ##     title(paste("Depth=", Apf$H$z))
dev.set( g$MAINdev)
  }
