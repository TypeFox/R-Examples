GEOmap.Extract<-function(MAP, SEL, INOUT="out")
{

    if(missing(INOUT)) INOUT ="out"

    INOUT =tolower(INOUT)
    
###    SEL 
  GM = GEOmap.list(MAP)

    
 useit = 1:length(GM$STROKES$nam)

  if(INOUT=="in")
    {
     useit =  SEL
    }
    else
      {

        useit =  useit[-SEL]

      }
    
  NEWMAP = list(STROKES = list(nam = NULL, num = NULL, index = NULL,
                  col = NULL, style = NULL, code = NULL, LAT1 = NULL, LAT2 = NULL,
                  LON1 = NULL, LON2 = NULL), POINTS = list(lat = NULL,
                                               lon = NULL), LL=list())
 

  kmap = 1
  for(i in useit)
    {
      NEWMAP$STROKES$nam[kmap] =  GM$STROKES$nam[i]
      NEWMAP$STROKES$num[kmap] =  GM$STROKES$num[i]
      NEWMAP$STROKES$index[kmap] = GM$STROKES$index[i]
      NEWMAP$STROKES$col[kmap] =  GM$STROKES$col[i] 
      NEWMAP$STROKES$style[kmap] = GM$STROKES$style[i]
      NEWMAP$STROKES$code[kmap] = GM$STROKES$code[i]
      NEWMAP$STROKES$LAT1[kmap] = GM$STROKES$LAT1[i]
      NEWMAP$STROKES$LAT2[kmap] = GM$STROKES$LAT2[i]
      NEWMAP$STROKES$LON1[kmap] = GM$STROKES$LON1[i]
      NEWMAP$STROKES$LON2[kmap] = GM$STROKES$LON2[i]
      
      NEWMAP$LL[[kmap]] =  GM$LL[[i]]
      kmap = kmap+1
    }

  GM = list.GEOmap(NEWMAP)

  
  invisible(GM)



}


