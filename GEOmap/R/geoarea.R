
geoarea<-function(MAP, proj=NULL , ncut=10)
  {
    ########  calculate the area of map elements in GEOmap
    if(missing(ncut)) ncut=10 
    if(missing(proj)) proj=NULL 

   ##  if(!require(splancs)) { print("NEED splancs"); return(NULL) }

    cgAREA=rep(NA,length(MAP$STROKES$num)) 

    for(i in 1:length(MAP$STROKES$num))
      {

        if(MAP$STROKES$num[i] > ncut)
          {
            j1 = MAP$STROKES$index[i] + 1
            j2 = j1 + MAP$STROKES$num[i] - 1
            JEC = j1:j2
            lon = MAP$POINTS$lon[JEC]
            lat = MAP$POINTS$lat[JEC]

            if(is.null(proj))
              {
                PROJ = setPROJ(type=2, LAT0=median(lat), LON0=median(lon))
              }
            else
              {
                PROJ = proj

              }

            X = GLOB.XY(lat, lon, PROJ)
            POL = cbind(X$x, X$y)

            cgAREA[i] = splancs::areapl(POL)

          }
        else
          {
            cgAREA[i] = NA
          }




      }

    return(cgAREA)


  }
