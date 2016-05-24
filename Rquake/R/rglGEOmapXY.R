`rglGEOmapXY` <-
function(MAP, LIM=c(-180, -90, 180, 90), PROJ=list(),  PMAT=NULL,
         add=TRUE,  SEL=NULL, GRID=NULL, GRIDcol=1,
         MAPcol=NULL, MAPstyle=NULL, border=NA, cenlon=0, shiftlon=0,
         linelty=1, linelwd=1, ptpch=".", ptcex=1, NUMB=FALSE, ...)
{
  ###  NUMB = add a number on the stroke to show which stroke it is for later modification
  ######   MAPcol will override the color in the map data base.  Good for making BW figures
  if(missing(cenlon)) { cenlon=0 }
  if(missing(GRID)) { GRID=NULL }
    if(missing(GRIDcol)) { GRIDcol=1 }
  if(missing(PMAT)) { PMAT=NULL }
  if(missing(SEL)) {  SEL=NULL }
  if(missing(linelty)) { linelty=1 }
  if(missing(NUMB)) { NUMB=FALSE }
  if(missing(linelwd)) { linelwd=1 }
  if(missing(ptpch)) { ptpch="." }
  if(missing(ptcex)) { ptcex=1 }

  
  if(missing(MAPcol)) { MAPcol=NULL }
  if(missing(MAPstyle)) { MAPstyle=NULL }
  if(missing(border)) { border=NA  }
  if(missing(shiftlon)) {  shiftlon=0 }

  
  pctexp = 0.01
 ###  MAP$POINTS$lon = RPMG::fmod( MAP$POINTS$lon, 360)

  ###  LIM can be a vector (lon1, lat1, lon2, lat2)
  ###  or LIM can be a list: list(lon=c(lon1, lon2), lat=c(lat1, lat2))

  
  if(missing(LIM))
    {
      lon = RPMG::fmod(MAP$POINTS$lon-shiftlon, 360)
      
     ###  LIMP=c( min(lon), min(MAP$POINTS$lat), max(lon), max(MAP$POINTS$lat))

      RLON = GEOmap::expandbound( range(RPMG::fmod(MAP$POINTS$lon, 360)), pctexp) 
      RLAT = GEOmap::expandbound( range(MAP$POINTS$lat), pctexp) 

      
      LIM=c( RLON[1], RLAT[1], RLON[2] , RLAT[2] )

    
    }
  else
    {
      if(is.null(LIM))
        {
          RLON = GEOmap::expandbound( range(RPMG::fmod(MAP$POINTS$lon, 360)), pctexp) 
          RLAT = GEOmap::expandbound( range(MAP$POINTS$lat), pctexp) 
          
          LIM=c( RLON[1], RLAT[1], RLON[2] , RLAT[2] )

        }

      
      if(is.list(LIM))
        {
         
          lon = RPMG::fmod(LIM$lon-shiftlon, 360)
          lat = LIM$lat
         ###  LIMP=c( min(lon), min(lat), max(lon), max(lat))
          LIM=c( min(RPMG::fmod(lon, 360)), min(lat), max(fmod(lon, 360)), max(lat))
          
        }
    }
  
  

  LLlim = list(lat=LIM[c(2,4)], lon=LIM[c(1,3)])

########################################################
  if(missing(PROJ)) {
    ############   need a projection  ###############
   ### MAPMEDLAT =  median(MAP$POINTS$lat)
   ### MAPMEDLON =  median(MAP$POINTS$lon-shiftlon)

    MAPCENLAT =  mean(LLlim$lat)
    MAPCENLON =  median(LLlim$lon)

    PROJ = GEOmap::setPROJ(type=2, LAT0=MAPCENLAT, LON0=MAPCENLON ,
      LATS=NULL, LONS=NULL, DLAT=NULL, DLON=NULL, FN =0)
    
  }




  
  
  if(missing(add)) { add=FALSE }

  if(is.null(MAP$POINTS$z))
     {

       MAP$POINTS$z = rep(0, length(MAP$POINTS$lat))

     }

  if(is.null(MAP$POINTS$x))
    {
      MAPXY = GEOmap::GLOB.XY(MAP$POINTS$lat ,  RPMG::fmod( MAP$POINTS$lon-shiftlon, 360) , PROJ )



      
      MAXDISTX = max( MAPXY$x) - min(MAPXY$x)
      
      if(is.null(PMAT))
        {
          MAP$POINTS$x = MAPXY$x
          MAP$POINTS$y = MAPXY$y
        }
      else
        {
          
          tem = trans3d(MAPXY$x, MAPXY$y, MAP$POINTS$z , PMAT)
          MAP$POINTS$x = tem$x
          MAP$POINTS$y = tem$y
        }

      STRKXYLL = GEOmap::GLOB.XY( MAP$STROKES$LAT1,  RPMG::fmod(MAP$STROKES$LON1-shiftlon, 360)  , PROJ )
      STRKXYUR = GEOmap::GLOB.XY( MAP$STROKES$LAT2,  RPMG::fmod(MAP$STROKES$LON2-shiftlon, 360)  , PROJ )

      STRKXYUL = GEOmap::GLOB.XY( MAP$STROKES$LAT2,  RPMG::fmod(MAP$STROKES$LON1-shiftlon, 360)  , PROJ )
      STRKXYLR = GEOmap::GLOB.XY( MAP$STROKES$LAT1,  RPMG::fmod(MAP$STROKES$LON2-shiftlon, 360)  , PROJ )


      
      if(is.null(PMAT))
        {
          MAP$STROKES$x1 = STRKXYLL$x
          MAP$STROKES$y1 = STRKXYLL$y
        }
      else
        {
          tem = trans3d(STRKXYLL$x, STRKXYLL$y, rep(0, length(STRKXYLL$y)) , PMAT)

          MAP$STROKES$x1= tem$x
          MAP$STROKES$y1= tem$y
        }


      
      if(is.null(PMAT))
        {
          MAP$STROKES$x2 = STRKXYUR$x
          MAP$STROKES$y2 = STRKXYUR$y
        }
      else
        {

          tem = trans3d(STRKXYUR$x, STRKXYUR$y, rep(0, length(STRKXYUR$y)) , PMAT)

          MAP$STROKES$x2= tem$x
          MAP$STROKES$y2= tem$y
        }

    }
  
 ##  print(paste(sep=" ", "test", min(MAP$POINTS$x)))
 ##  print(LIM)


      XYLIM =  GEOmap::GLOB.XY(LLlim$lat,LLlim$lon, PROJ)
      LLlim$x = XYLIM$x
      LLlim$y = XYLIM$y

  
  xrange = abs( diff(range(LLlim$x, na.rm=TRUE)) ) 
  yrange = abs( diff(range(LLlim$y, na.rm=TRUE)) )

 ## print(c(xrange, yrange))
  

  Kstroke = length(MAP$STROKES$num)


  FORCE = FALSE

  if(!is.null(SEL))
    {
      if(SEL<0)
        {
          FORCE = TRUE
          SEL = 1:length(MAP$STROKES$num)
        }
    }

  ###  if(exists("worldmap"))

  if(FORCE==FALSE)
    {
      IN = GEOmap::KINOUT(MAP, LLlim , projtype=2)
    }
  else
    {
      ###########  plot all the strokes  
      IN = 1:length(MAP$STROKES$num)

    }


##########    for selecting parts of the map
  if(!is.null(SEL))
    {

      slin = which(IN %in% SEL)
      if(length(slin)>=1)
        { IN = IN[slin] }
    }
 

  

  minx=Inf; maxx=-Inf; miny=Inf; maxy=-Inf;
  
  ##  print(IN)
  if(length(IN)<1)
    {
      print("No map strokes in target")
      return(0)

    }

##############   this inserts NA where the lines cross the boundaries 
  for(i in IN)
    {
      
      j1 = MAP$STROKES$index[i]+1
      j2 = j1+MAP$STROKES$num[i]-1
      
      if(j1>0 & j2>0 & j2-j1 >=0)
        {
          JEC = j1:j2
          x = MAP$POINTS$x[JEC]
          y = MAP$POINTS$y[JEC]
          x[x<LLlim$x[1] |  x>LLlim$x[2] ] = NA
          y[y<LLlim$y[1] |  y>LLlim$y[2] ] = NA
          minx = min(c(minx,x), na.rm=TRUE)
          maxx = max(c(maxx,x), na.rm=TRUE)
          miny = min(c(miny,y), na.rm=TRUE)
          maxy = max(c(maxy,y), na.rm=TRUE)
        }
      else
        {
          next
        }
      
    }

##############   this replaces the colors with a fixed color
  if(!is.null(MAPcol))
    {

      MAP$STROKES$col = rep(MAPcol, length=length(MAP$STROKES$col))

    }
  if(!is.null(MAPstyle))
    {

      MAP$STROKES$style = rep(MAPstyle, length=length(MAP$STROKES$style))

    }

  ##############

  if(add==FALSE)
    {
###R1 = range(c(MAP$STROKES$x1[IN],MAP$STROKES$x2[IN]))
###R2 = range(c(MAP$STROKES$y1[IN], MAP$STROKES$y2[IN]))
 ###     plot(c(minx, maxx) , c(miny,maxy), asp=1, type='n', ...)

 plot(range(LLlim$x), range(LLlim$y),  asp=1, type='n', ...)
      
      
### plot(c(MAP$STROKES$x1[IN],MAP$STROKES$x2[IN]) , c(MAP$STROKES$y1[IN], MAP$STROKES$y2[IN]), asp=TRUE, type='n', ...)
###print(c(R1, R2))
    }

  
  if(!is.null(GRID))
    {
      GEOmap::addLLXY(GRID$lats, GRID$lons, PMAT=PMAT, GRIDcol=GRIDcol, LABS=0, BORDER=0 , PROJ=PROJ )
    }
  
  for(i in IN)
    {
      j1 = MAP$STROKES$index[i]+1
      j2 = j1+MAP$STROKES$num[i]-1



 ###cat(paste(i, "style=",  MAP$STROKES$style[i]), sep="\n")
      
      if( (j1>0 & j2>0 & j2-j1 >=0))
        {
          JEC = j1:j2
        }
      else
        {
          next
        }
      
      ###  print(paste(sep=' ',"----------------------", i, j1,j2))

      if(NUMB==TRUE)
        {
         ### points(MAP$POINTS$x[j1], MAP$POINTS$y[j1])
         ### text(MAP$POINTS$x[j1], MAP$POINTS$y[j1], labels=i, pos=3)
        }

      
      if(MAP$STROKES$style[i]==1)
        {
         ### cat(paste(i, "***********plotting style=1"), sep="\n")

          zp =rep(0, length=length(MAP$POINTS$y[JEC]))
          
          rgl::points3d(cbind(MAP$POINTS$x[JEC], MAP$POINTS$y[JEC],zp), col=MAP$STROKES$col[i] )

          
        }

      
      if(MAP$STROKES$style[i]==2)
        {
          x = MAP$POINTS$x[JEC]
          y = MAP$POINTS$y[JEC]

         
          ###############################  this code is meant to
          ###############################  to avoid wraparound
          xd = abs(diff(c(x, x[1])))
          wwx = which(xd>0.9*xrange)

          yd = abs(diff(c(y, y[1])))
          wwy = which(yd>0.80*yrange)

          ###  print(wwx) 
          ###  print(wwy) 

          
          if( (!is.null(wwx) & length(wwx)>0) | (!is.null(wwy) & length(wwy)>0) )
            {
             ## print(paste(sep=' ', "################", i, length(x), length(ww)))
              
              ##  print(ww)
              if(length(wwx)>0)
                {
                  zy = GEOmap::insertNA(y, wwx)
                  zx = GEOmap::insertNA(x, wwx)
                  zz = rep(0, length(zx))
                 rgl::lines3d(cbind(zx, zy, zz) , col=MAP$STROKES$col[i] )

                 ### lines(zx, zy, col='blue' , lty=linelty, lwd=linelwd)
            

                }


              if(length(wwy)>0)
                {
                  zy = GEOmap::insertNA(y, wwy)
                  zx = GEOmap::insertNA(x, wwy)
                  zz = rep(0, length(zx))
                 rgl::lines3d(cbind(zx, zy, zz) , col=MAP$STROKES$col[i] )

                  

                 ###  lines(zx, zy, col='blue' , lty=linelty, lwd=linelwd)
                  
                }



              

              
            }
          else
            {
             ## x[MAP$POINTS$lon[JEC]<LLlim$lon[1] |  MAP$POINTS$lon[JEC]>LLlim$lon[2] ] = NA
             ## y[MAP$POINTS$lat[JEC]<LLlim$lat[1] |  MAP$POINTS$lat[JEC]>LLlim$lat[2] ] = NA

            ##  print(LLlim)
              
             ## print(cbind(x,y)) 
              
              x[MAP$POINTS$x[JEC]<LLlim$x[1] |  MAP$POINTS$x[JEC]>LLlim$x[2] ] = NA
              y[MAP$POINTS$y[JEC]<LLlim$y[1] |  MAP$POINTS$y[JEC]>LLlim$y[2] ] = NA


              
              ##   lines(x, y, col='blue' )

            ##  print(cbind(x,y)) 
               zz = rep(0, length(x))
                 rgl::lines3d(cbind(x, y, zz) , col=MAP$STROKES$col[i] )

             ##   lines(x, y, col=MAP$STROKES$col[i], lty=linelty, lwd=linelwd)
            } 
        }

      if(MAP$STROKES$style[i]==3)
        {

          x = MAP$POINTS$x[JEC]
          y = MAP$POINTS$y[JEC]


          x = c(x, x[1])
          y = c(y, y[1])

          x[MAP$POINTS$lon[JEC]<LLlim$lon[1] |  MAP$POINTS$lon[JEC]>LLlim$lon[2] ] = NA
          y[MAP$POINTS$lat[JEC]<LLlim$lat[1] |  MAP$POINTS$lat[JEC]>LLlim$lat[2] ] = NA

           ##    polygon(x, y, border=border, col=MAP$STROKES$col[i])
          
        ##   MM = fixCoastwrap(list(x=x, y=y), MAXDISTX )
        ##   polygon(MM$x, MM$y, border=border, col=MAP$STROKES$col[i])



          
        }



      
      ##  locator()

    }

  invisible(IN)
  

}

