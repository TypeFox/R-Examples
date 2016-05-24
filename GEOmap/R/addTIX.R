`addTIX` <-
function(lats, lons,   PROJ=list(), PMAT=NULL, col=gray(0.7), TICS=c(1,1), OUTER=TRUE, sides=c(1,2,3,4))
{
   if(missing(col)) {col=gray(0.7)}
   if(missing(TICS)) { TICS=NULL }
   if(missing(OUTER)) { OUTER=TRUE }
   if(missing(sides)) { sides=c(1,2,3,4)  }
   if(missing(PMAT)) { PMAT=NULL  }


mypar = par("usr")
   

   #####   first do longitudes along a given latitude:
   if(is.null(TICS)) { return() }

   #####

   if(any(sides==1))
     {
        i = lats[1]

   if(OUTER==TRUE)
     {
        j = lons
      }
   else
     {
       j = lons[-c(1, length(lons))]
     }
        
        xy = GLOB.XY(rep(i, length(j)), j ,  PROJ)
        txy = GLOB.XY(rep(i, length(j))+TICS[1], j , PROJ )
        if(!is.null(PMAT))
          {
            tem = trans3d(xy$x, xy$y, rep(0, length(xy$y)) , PMAT)
            xy$x = tem$x
            xy$y = tem$y
            tem = trans3d(txy$x, txy$y, rep(0, length(txy$y)) , PMAT)
            txy$x = tem$x
            txy$y = tem$y
          }

        segments(xy$x, xy$y, txy$x, txy$y)
        ## print(xy$x)
        ##  print(xy$y)
      }

   if(any(sides==3))
     {
       i = lats[length(lats)]
       
       if(OUTER==TRUE)
         {
           j = lons
         }
       else
         {
           j = lons[-c(1, length(lons))]
         }
       
       xy = GLOB.XY(rep(i, length(j)), j , PROJ )
       txy = GLOB.XY(rep(i, length(j))-TICS[1], j  , PROJ)
       if(!is.null(PMAT))
         {
           tem = trans3d(xy$x, xy$y, rep(0, length(xy$y)) , PMAT)
           xy$x = tem$x
           xy$y = tem$y
           tem = trans3d(txy$x, txy$y, rep(0, length(txy$y)) , PMAT)
            txy$x = tem$x
           txy$y = tem$y
         }
       
       segments(xy$x, xy$y, txy$x, txy$y)
       
     }


   ###  then do latitudes along a given longitude
   if(any(sides==2))
     {
       if(OUTER==TRUE)
         {        
           i = lats
         }
       else
         {
           i = lats[-c(1, length(lats))]
           
         }
       j = lons[1]
       
       xy = GLOB.XY(i,  rep(j, length(i)), PROJ)
       txy = GLOB.XY(i,  rep(j, length(i))+TICS[2], PROJ )
       
       if(!is.null(PMAT))
         {
           tem = trans3d(xy$x, xy$y, rep(0, length(xy$y)) , PMAT)
           xy$x = tem$x
           xy$y = tem$y
           
           tem = trans3d(txy$x, txy$y, rep(0, length(txy$y)) , PMAT)
           txy$x = tem$x
           txy$y = tem$y
           
         }
       segments(xy$x, xy$y, txy$x, txy$y)
     }
   if(any(sides==4))
     {
       if(OUTER==TRUE)
         {        
           i = lats
         }
       else
         {
           i = lats[-c(1, length(lats))]
           
         }
       j = lons[length(lons)]
       
       xy = GLOB.XY(i,  rep(j, length(i)) ,PROJ )
       txy = GLOB.XY(i,  rep(j, length(i))-TICS[2], PROJ )
       
       if(!is.null(PMAT))
         {
           tem = trans3d(xy$x, xy$y, rep(0, length(xy$y)) , PMAT)
           xy$x = tem$x
           xy$y = tem$y
           
           tem = trans3d(txy$x, txy$y, rep(0, length(txy$y)) , PMAT)
           txy$x = tem$x
           txy$y = tem$y
           
         }
       segments(xy$x, xy$y, txy$x, txy$y)
     }
      
}

