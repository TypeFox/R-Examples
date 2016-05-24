mapTeleSeis <-
function(sta, mylist, worldmap=NULL)
  {
   
    stinfo = list(mlat=median(sta$lat), mlon=median(sta$lon) ) 
    
    plotGEOmap(worldmap ,MAPcol=grey(.75),  add=FALSE, xaxs='i')

    points(  RPMG::fmod(sta$lon, 360)    ,  sta$lat, col='red', pch=6)

    pcol = rainbow(length(mylist))
    
    for(i in 1:length(mylist))
      {
        gev = mylist[[i]]
        points(  RPMG::fmod(gev$lon, 360)    ,  gev$lat, col='purple', pch=8)
        text( RPMG::fmod(gev$lon, 360)    ,  gev$lat, labels=i, pos=3)
        garc = getgreatarc( gev$lat  ,RPMG::fmod(gev$lon, 360) ,  stinfo$mlat , fmod(stinfo$mlon,360) , num = 50)

        flon = RPMG::fmod(garc$lon, 360)
        flat = garc$lat

##  print(i)
      ##  print(flon)
        
        if(any(abs(diff(flon))>100))
          {
            wsplit = which(abs(diff(flon))>100)
           #   indde = 1:length(flon)
            flon =  insertNA(flon,wsplit )
            flat =  insertNA(flat,wsplit )

           
           # indde =  insertNA(indde,wsplit )

           ## text(flon, flat, labels=indde, pos=1)

           ## print(cbind(flon, flat, indde))
          }


     ##    print(flon)
     ##    print(flat)

        
     #   points(flon    ,  flat , col= pcol[i] )
  lines(flon    ,  flat , col= pcol[i] )
        
    ##    lines( Fline$x     ,  Fline$y  , col='green' )
      }
  }

