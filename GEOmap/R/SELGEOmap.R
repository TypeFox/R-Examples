SELGEOmap<-function (MAP, ncut = 3, acut = c(0, 1e+05), proj = NULL, LIM = NULL) 
{
  if (missing(ncut)) 
    ncut = 0
  if (missing(proj)) 
    proj = NULL
  if (missing(acut)) 
    acut = NULL
  if (missing(LIM)) 
    LIM = NULL
 # if (!require(splancs)) {
 #   print("NEED splancs")
#    return(NULL)
#  }
  NEWMAP = list(STROKES = list(nam = NULL, num = NULL, index = NULL, 
                  col = NULL, style = NULL, code = NULL, LAT1 = NULL, LAT2 = NULL, 
                  LON1 = NULL, LON2 = NULL), POINTS = list(lat = NULL, 
                                               lon = NULL))
  
  if(is.null(proj))
    {
      proj = setPROJ(type = 1, LAT0 = 0, LON0 = 0, LAT1 = 0, LAT2 = 0, LATS = NULL, LONS = NULL, DLAT = NULL, DLON = NULL, FE = 0, FN = 0)
    }

  
   if (is.null(LIM))
     {

       cat("NO GEOGRAPHIC limits provided:", sep="\n"  )
       cat("Using Default:  LIM = c(0.0000, -70 ,  359.9167,   80 )   ", sep="\n"  )
       
       LIM = c(0.0000, -70 ,  359.9167,   80 ) 
       
       LLlim = list(lat = sort(LIM[c(2, 4)]), lon = sort(LIM[c(1, 3)]) )
       proj = setPROJ(type = 1, LAT0 = 0, LON0 = 0, LAT1 = 0, LAT2 = 0, LATS = NULL, LONS = NULL, DLAT = NULL, DLON = NULL, FE = 0, FN = 0)
      
       
       
     }
  else
    {

      LLlim = list(lat = sort(LIM[c(2, 4)]), lon = sort(LIM[c(1, 3)]) )
     
      
    }

    LIMXY = GLOB.XY(LLlim$lat , LLlim$lon, proj)
  
  
  Kmap = 0
  index1 = 0
  
  
  useit = FALSE
  OKarea = FALSE
  OKnum = FALSE


  for(i in 1:length(MAP$STROKES$num))
    {
      useit = FALSE
      OKarea = FALSE
      OKnum = FALSE
      INside = FALSE
      if (is.null(acut)) {OKarea = TRUE }
      if (is.null(ncut)) {OKnum = TRUE }


      
      if (MAP$STROKES$num[i] > ncut)
        {
          j1 = MAP$STROKES$index[i] + 1
          j2 = j1 + MAP$STROKES$num[i] - 1
          JEC = j1:j2
          lon = MAP$POINTS$lon[JEC]
          lat = MAP$POINTS$lat[JEC]
          

          PROJ = setPROJ(type = 2, LAT0 = median(lat), 
            LON0 = median(lon))
          
          X = GLOB.XY(lat, lon, PROJ)
          POL = cbind(X$x, X$y)
          
          OKnum = TRUE
          if (!is.null(acut)) {
            cgAREA = splancs::areapl(POL)
            if (cgAREA >= acut[1] & cgAREA <= acut[2]) {
              OKarea = TRUE
            }
          }

          if(OKarea &  OKnum)
            {

              LL2XY = GLOB.XY(lat, lon, proj)

              if(any(LL2XY$x>=LIMXY$x[1] & LL2XY$x<=LIMXY$x[2]
                     & LL2XY$y>=LIMXY$y[1] & LL2XY$y<=LIMXY$y[2] ))
                {
                  INside = TRUE
                }
              else
                {
                  INside = FALSE
                }

            }
          
        }
      useit = OKarea & OKnum & INside
      if(useit) {
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
