GEOmap.cat<-function(MAP1, MAP2)
{

  ## Combine two GEOmap list structures

  G1 = GEOmap.list(MAP1)

  G2 = GEOmap.list(MAP2)

  G3 = list(STROKES = list(nam = NULL, num = NULL, index = NULL,
              col = NULL, style = NULL, code = NULL, LAT1 = NULL, LAT2 = NULL,
              LON1 = NULL, LON2 = NULL), POINTS = list(lat = NULL,
                                           lon = NULL), LL=list() )
  
################   G3$STROKES = c(G1$STROKES$, G2$STROKES$)


  G3$STROKES =
    list(nam =c(G1$STROKES$nam, G2$STROKES$nam) ,
         num =c(G1$STROKES$num, G2$STROKES$num) ,
         index =c(G1$STROKES$index, G2$STROKES$index) ,
         col =c(G1$STROKES$col, G2$STROKES$col) ,
         style =c(G1$STROKES$style, G2$STROKES$style) ,
         code =c(G1$STROKES$code, G2$STROKES$code) ,
         LAT1 =c(G1$STROKES$LAT1, G2$STROKES$LAT1) ,
         LAT2 =c(G1$STROKES$LAT2, G2$STROKES$LAT2) ,
         LON1 =c(G1$STROKES$LON1, G2$STROKES$LON1) ,
         LON2 =c(G1$STROKES$LON2, G2$STROKES$LON2) )


  G3$LL = c(G1$LL, G2$LL)

###   print(c(length(G1$LL), length(G2$LL), length(G3$LL) ))

  GOUT = list.GEOmap(G3) 

  invisible(GOUT) 

}
