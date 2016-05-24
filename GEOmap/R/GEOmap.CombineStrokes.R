GEOmap.CombineStrokes<-function(MAP, SEL)
{
###    SEL should be in order of the combination
  GM = GEOmap.list(MAP)

  j = SEL[1]

  newstroke = list(x=NULL, y=NULL)
  
  for(i in 1:length(SEL))
    {
       newstroke$x = c(newstroke$x, GM$LL[[ SEL[i] ]]$lon)
       newstroke$y = c(newstroke$y,  GM$LL[[ SEL[i] ]]$lat)
    }


  GM$LL[[ j ]]$lon = newstroke$x
  GM$LL[[ j ]]$lat = newstroke$y
  

  NEWMAP = list(STROKES = list(nam = NULL, num = NULL, index = NULL,
                  col = NULL, style = NULL, code = NULL, LAT1 = NULL, LAT2 = NULL,
                  LON1 = NULL, LON2 = NULL), POINTS = list(lat = NULL,
                                               lon = NULL))


  rem = SEL[2:length(SEL)]


  CHOOS =  1:length(MAP$STROKES$nam)
  CHOOS = CHOOS[-rem]

  kmap = 1
  for(i in CHOOS)
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


