DiagnoseFailure <- structure(function# Visualize points that could not be mapped using RapidPolygonLookup()
                             ### This functions plots the points that could not be mapped using RapidPolygonLookup()
                             ### The points are overlayed on the polygons to contextualize their
                             ### geographical location and understand the reason behind their exclusion.
                             (XY.polys, ##<< output from function RapidPolygonLookup()
                              poly.list= NULL ##<< polygon list with 3 or 4 elements: data, polys, poly.centers, and possibly ranges. Needs to be supplied if RapidPolygonLookup() was run with keep.data= FALSE
                              ){
                               still.unmapped <- which(is.na(XY.polys$XY[,3]))
                               if (length(still.unmapped) == 0){
                                 return(cat("No unmapped points found in data \n"))
                               }
                               # diagnostic plotting code to visually understand why these points are so difficult to map.
                               XLIM <- range(XY.polys$XY[still.unmapped, c("X")])
                               YLIM <- range(XY.polys$XY[still.unmapped, c("Y")])
                               #expand the plot range by 100%:
                               XLIM <- mean(XLIM) + diff(XLIM)*c(-1, 1)*2
                               YLIM <- mean(YLIM) + diff(YLIM)*c(-1, 1)*2
                               plotPolys(XY.polys$poly.list$polys, xlim= XLIM, ylim= YLIM)
                               points(XY.polys$XY[still.unmapped,c("X","Y")], col="red")
                               points(XY.polys$poly.list$poly.centers[,c("X","Y")], col="green")
                               
                               for (i in still.unmapped) {
                                 points(XY.polys$poly.list$poly.centers, col="green")
                               }
                               # Plot of the points left unmapped after RapidPolygonLookup call
                             }, ex= function(){
                               data(sf.crime.2012, envir = environment())
                               data(sf.polys, envir = environment())
                               cat(nrow(sf.crime.2012), "rows in SF crime \n")
                               
                               XY.kdtree <- RapidPolygonLookup(sf.crime.2012[,c("X","Y")], poly.list= sf.polys, 
                                                               k= 10, N= 1000, 
                                                               poly.id= "fips", poly.id.colname= "census.block", 
                                                               keep.data= TRUE, verbose= TRUE)
                               DiagnoseFailure(XY.kdtree)
                             }
  )
  
  