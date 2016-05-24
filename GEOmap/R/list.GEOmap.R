list.GEOmap<-function(MAP, SEL=1)
  {
########   return a GEOmap list to a normal GEOmap with points as one long vector
    ########

    N = length(MAP$LL)
    
    if(missing(SEL)) SEL = 1:N 
    
    NEWMAP = list(STROKES = list(nam = NULL, num = NULL, index = NULL,
                    col = NULL, style = NULL, code = NULL, LAT1 = NULL, LAT2 = NULL,
                    LON1 = NULL, LON2 = NULL), POINTS = list(lat = NULL,
                                                 lon = NULL))
    Kmap = 1
    index1 = 0

    lon = NULL
    lat = NULL
NEWMAP$POINTS$lon=NULL
NEWMAP$POINTS$lat=NULL

    if(is.null(MAP$STROKES))
      {
        LON2=LON1=LAT2=LAT1=vector(); 
        for(i in 1:N)
          {
            LAT1[i] = min(MAP$LL[[i]]$lat)
               LAT2[i]  = max(MAP$LL[[i]]$lat)
               LON1[i]  = min(MAP$LL[[i]]$lon)
               LON2[i]  = max(MAP$LL[[i]]$lon)


          }

        
        MAP$STROKES = 
          list(nam =paste("m",1:N, sep="")  ,
               num =unlist(lapply(MAP$LL, length)) ,
               index = NULL,
               col = rep("black", N),
               style = rep(2, N),
               code = rep("o", N),
               LAT1 = LAT1,
               LAT2 = LAT2,
               LON1 = LON1,
               LON2 = LON2)
        
      }
    
    for(i in SEL)
      {

        lon = MAP$LL[[i]]$lon
        lat = MAP$LL[[i]]$lat
        
        NEWMAP$STROKES$nam[Kmap] =  MAP$STROKES$nam[i]
        NEWMAP$STROKES$num[Kmap] = length(lon)
        NEWMAP$STROKES$index[Kmap] = index1
        NEWMAP$STROKES$col[Kmap] =  MAP$STROKES$col[i] 
        NEWMAP$STROKES$style[Kmap] = MAP$STROKES$style[i] 
        NEWMAP$STROKES$code[Kmap] = MAP$STROKES$code[i] 
        NEWMAP$STROKES$LAT1[Kmap] = MAP$STROKES$LAT1[i] 
        NEWMAP$STROKES$LAT2[Kmap] = MAP$STROKES$LAT2[i] 
        NEWMAP$STROKES$LON1[Kmap] = MAP$STROKES$LON1[i] 
        NEWMAP$STROKES$LON2[Kmap] = MAP$STROKES$LON2[i] 
        NEWMAP$POINTS$lon = c(NEWMAP$POINTS$lon,  lon)
        NEWMAP$POINTS$lat = c(NEWMAP$POINTS$lat,   lat)
        index1 = index1+length(lon)
        Kmap = Kmap + 1

        
      }

    return(NEWMAP)
  }


