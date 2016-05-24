GEOmap.limit<-function(MAP, LLlim )
{

  ###   LLlim =c(lon1, Lat1, Lon2, Lat2)

  if(is.list(LLlim))
    {
      lat1 = LLlim$lat[1]
      lon1= RPMG::fmod(LLlim$lon[1], 360)
      lat2= LLlim$lat[2]
      lon2= RPMG::fmod(LLlim$lon[2], 360)

    }
  else
    {
      lat1 = LLlim[2]
      lon1= RPMG::fmod(LLlim[1], 360)
      lat2= LLlim[4]
      lon2= RPMG::fmod(LLlim[3], 360)
    }

  
###    SEL 
  GM = GEOmap.list(MAP)

    
 useit = 1:length(GM$STROKES$nam)

  NEWMAP = list(STROKES = list(nam = NULL, num = NULL, index = NULL,
                  col = NULL, style = NULL, code = NULL, LAT1 = NULL, LAT2 = NULL,
                  LON1 = NULL, LON2 = NULL), POINTS = list(lat = NULL,
                                               lon = NULL), LL=list())
 

  kmap = 0
  index1 = 0

  for(i in useit)
    {
      LL =  GM$LL[[i]]

      GetTF = LL$lat>lat1 & LL$lon>lon1 & LL$lat<lat2 & LL$lon<lon2

      WLL = which(GetTF)

      ##  GEOmap.breakline(Y, 5)
      if(length(WLL)<1) next

      dWLL = diff(WLL)
      nix = which(dWLL>1)
      kix = length(nix)
      if(kix>0)
        {
        ##  print(c(i, kix, nix))
          Zbreak  = RPMG::breakline.index(WLL, nix)
         ## print(length(Zbreak))
          
          for(ib in 1:length(Zbreak))
            {
              kmap = kmap+1
              indez= Zbreak[[ib]]
              lat = LL$lat[indez]
              lon = LL$lon[indez]
              Npts = length(lon )
              RElatlon = list(lat=lat, lon=lon )
              
              NEWMAP$STROKES$nam[kmap] =  GM$STROKES$nam[i]
              NEWMAP$STROKES$num[kmap] =  Npts 
              index1 = index1+Npts 
              
              NEWMAP$STROKES$index[kmap] = index1
              NEWMAP$STROKES$col[kmap] =  GM$STROKES$col[i] 
              ## NEWMAP$STROKES$style[kmap] = GM$STROKES$style[i]
                                        # for now, change the style to a stroke
              NEWMAP$STROKES$style[kmap] = 2
              NEWMAP$STROKES$code[kmap] = GM$STROKES$code[i]
              NEWMAP$STROKES$LAT1[kmap] = min(RElatlon$lat)
              NEWMAP$STROKES$LAT2[kmap] = max(RElatlon$lat)
              NEWMAP$STROKES$LON1[kmap] = min(RElatlon$lon)
              NEWMAP$STROKES$LON2[kmap] = max(RElatlon$lon)
           
              NEWMAP$LL[[kmap]] =  RElatlon
            }
        }
      else
        {
      kmap = kmap+1
      Npts = length(LL$lon[WLL] )
      RElatlon = list(lat=LL$lat[WLL], lon=LL$lon[WLL] )
      
      NEWMAP$STROKES$nam[kmap] =  GM$STROKES$nam[i]
      NEWMAP$STROKES$num[kmap] =  Npts 
      index1 = index1+Npts 
      
      NEWMAP$STROKES$index[kmap] = index1
      NEWMAP$STROKES$col[kmap] =  GM$STROKES$col[i] 
      ## NEWMAP$STROKES$style[kmap] = GM$STROKES$style[i]
                                        # for now, change the style to a stroke
      NEWMAP$STROKES$style[kmap] = 2
      NEWMAP$STROKES$code[kmap] = GM$STROKES$code[i]
      NEWMAP$STROKES$LAT1[kmap] = min(RElatlon$lat)
      NEWMAP$STROKES$LAT2[kmap] = max(RElatlon$lat)
      NEWMAP$STROKES$LON1[kmap] = min(RElatlon$lon)
      NEWMAP$STROKES$LON2[kmap] = max(RElatlon$lon)
      
      NEWMAP$LL[[kmap]] =  RElatlon
    }
    }

  GM = list.GEOmap(NEWMAP)

  
  invisible(GM)



}


