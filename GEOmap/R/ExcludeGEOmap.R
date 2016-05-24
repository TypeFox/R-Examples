ExcludeGEOmap<-function(MAP, SEL, INOUT="out")
{

  if(missing(INOUT)) INOUT ="out" 
  NEWMAP = list(STROKES = list(nam = NULL, num = NULL, index = NULL,
                  col = NULL, style = NULL, code = NULL, LAT1 = NULL, LAT2 = NULL,
                  LON1 = NULL, LON2 = NULL), POINTS = list(lat = NULL,
                                               lon = NULL), LL=list())

  Kmap = 0
  index1 = 0
  for(i in 1:length(MAP$STROKES$nam))
    {

      j1 = MAP$STROKES$index[i] + 1
      j2 = j1 + MAP$STROKES$num[i] - 1
      JEC = j1:j2
      lon = MAP$POINTS$lon[JEC]
      lat = MAP$POINTS$lat[JEC]
          

      if(INOUT=="in")
        {
           useit = any(!is.na(match(i, SEL)))
          
        }
      else
        {
          useit = all(is.na(match(i, SEL)))
          ##  useit = !useit
        }
      
      
  if (useit) {

    Kmap = Kmap + 1
    
    NEWMAP$STROKES$nam = c(NEWMAP$STROKES$nam, MAP$STROKES$nam[i])
    NEWMAP$STROKES$num = c(NEWMAP$STROKES$num, MAP$STROKES$num[i])
    NEWMAP$STROKES$index = c(NEWMAP$STROKES$index, index1)
    index1 = index1 + MAP$STROKES$num[i]
    NEWMAP$STROKES$col = c(NEWMAP$STROKES$col, MAP$STROKES$col[i])
    NEWMAP$STROKES$style = c(NEWMAP$STROKES$style, MAP$STROKES$style[i])
    NEWMAP$STROKES$code = c(NEWMAP$STROKES$code, MAP$STROKES$code[i])
    NEWMAP$STROKES$LAT1 = c(NEWMAP$STROKES$LAT1, MAP$STROKES$LAT1[i])
    NEWMAP$STROKES$LAT2 = c(NEWMAP$STROKES$LAT2, MAP$STROKES$LAT2[i])
    NEWMAP$STROKES$LON1 = c(NEWMAP$STROKES$LON1, MAP$STROKES$LON1[i])
    NEWMAP$STROKES$LON2 = c(NEWMAP$STROKES$LON2, MAP$STROKES$LON2[i])
    NEWMAP$POINTS$lon = c(NEWMAP$POINTS$lon, lon)
    NEWMAP$POINTS$lat = c(NEWMAP$POINTS$lat, lat)
  }

}



  
  invisible(NEWMAP)

}
