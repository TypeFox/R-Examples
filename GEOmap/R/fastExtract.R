fastExtract<-function(MAP, SEL, INOUT="out")
{

    if(missing(INOUT)) INOUT ="out"

    INOUT =tolower(INOUT)
    
###    SEL 
###  GM = GEOmap.list(MAP)

    
 useit = 1:length(MAP$STROKES$nam)

  if(INOUT=="in")
    {
     useit =  SEL
    }
    else
      {

        useit =  useit[-SEL]

      }
    
  NEWMAP = list(STROKES = list(
                  nam = NULL,
                  num = NULL,
                  index = NULL,
                  col = NULL,
                  style = NULL,
                  code = NULL,
                  LAT1 = NULL,
                  LAT2 = NULL,
                  LON1 = NULL,
                  LON2 = NULL),
    POINTS = list(lat = NULL, lon = NULL), LL=list())
 

  kmap = 1
  for(i in useit)
    {
      j1 = MAP$STROKES$index[i] + 1
      j2 = j1 + MAP$STROKES$num[i] - 1
      JEC = j1:j2
      lon = MAP$POINTS$lon[JEC]
      lat = MAP$POINTS$lat[JEC]
      
      if(j1>0 & j2>0 & j2-j1 >=0)
        {
          
          NEWMAP$STROKES$nam[kmap] =  MAP$STROKES$nam[i]
          NEWMAP$STROKES$num[kmap] =  MAP$STROKES$num[i]
          NEWMAP$STROKES$index[kmap] = MAP$STROKES$index[i]
          NEWMAP$STROKES$col[kmap] =  MAP$STROKES$col[i] 
          NEWMAP$STROKES$style[kmap] = MAP$STROKES$style[i]
          NEWMAP$STROKES$code[kmap] = MAP$STROKES$code[i]
          NEWMAP$STROKES$LAT1[kmap] = MAP$STROKES$LAT1[i]
          NEWMAP$STROKES$LAT2[kmap] = MAP$STROKES$LAT2[i]
          NEWMAP$STROKES$LON1[kmap] = MAP$STROKES$LON1[i]
          NEWMAP$STROKES$LON2[kmap] = MAP$STROKES$LON2[i]
          
      NEWMAP$LL[[kmap]] =  list(lat=lat, lon=lon)
          kmap = kmap+1
          
        }
      

      
    }
    
    GM = list.GEOmap(NEWMAP)
    
    
    invisible(GM)



  }


