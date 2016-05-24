RapidPolygonLookup <- structure(function# Efficient spatial polygon search using kd-trees.
                                ### Given spatial partitions such as census blocks, ZIP codes or police district boundaries, we are
                                ### frequently faced with the need to spatially aggregate data. 
                                ### Unless efficient data structures are used, this can be a daunting task. 
                                ### The operation point.in.polygon() from the package sp is computationally expensive.
                                ### Here, we exploit kd-trees as efficient nearest neighbor search algorithm 
                                ### to dramatically reduce the effective number of polygons being searched.
                                ### Points that are left unmapped are put through a linear search to find the
                                ### associated polygon.
                                (XY, ##<< data frame containing X-Y or (lon-lat, long-lat, longitude-latitude) columns
                                 polygons, ##<< polygons to crop and add poly centres
                                 poly.list= NULL, ##<< polygon list with three elements: data, polys, and poly.centers as output from function CropSpatialPolygonsDataFrame()
                                 k= 10, ##<< maximum number of near neighbours to compute. The default value is set to 10
                                 N= nrow(XY), ##<< number of rows of XY to search
                                 poly.id= "fips", ##<< column name in poly.list$data containing the polygon identifier
                                 poly.id.colname= "census.block", ##<< desired column name in the output data frame containing the polygon identifier
                                 keep.data= TRUE, ##<< retain polygon list and centers for future referece
                                 verbose=0 ##<< level of verbosity
                                ){
                                  
                                  if(length(grep("^X$", colnames(XY))) == 0 || 
                                       length(grep("^Y$", colnames(XY))) == 0){
                                    x <- which(grepl("lon|long|longitude", colnames(XY), ignore.case= TRUE))
                                    y <- which(grepl("lat|latitude", colnames(XY), ignore.case= TRUE))
                                    if((length(x) | length(y)) == 0){
                                      stop("X-Y columns not found in data")
                                    } else{
                                      colnames(XY)[c(x,y)] <- c("X", "Y") 
                                    }
                                  }
                                  
                                  XY <- XY[1:N, c("X", "Y")]
                                  
                                  if (is.null(poly.list)){
                                    bb <- qbbox(XY$Y, XY$X, 
                                                margin = list(m = c(1, 1, 1, 1), TYPE = c("perc")))
                                    bb <- do.call("data.frame", bb)
                                    names(bb) <- c("Y", "X")
                                    poly.list <- CropSpatialPolygonsDataFrame(x= polygons, bb= bb, verbose= 0)
                                  }
                                  
                                  k <- min(k, nrow(poly.list$poly.centers))
                                  XY <- SearchForPolygon(poly.list, XY, k, poly.id, poly.id.colname, verbose)
                                  if (any(is.na(XY[, poly.id.colname]))){
                                    if (verbose) cat("Mapping points missed in initial run using range-search \n")
                                    XY.na <- XY[which(is.na(XY[, poly.id.colname])),]
                                    XY.na <- FindPolygonInRanges(poly.list, XY.na, poly.id, poly.id.colname, verbose)
                                    XY <- rbind(XY[which(!(is.na(XY[, poly.id.colname]))),], XY.na)
                                    XY <- XY[order(XY$sort_col),]
                                    XY$sort_col <- NULL
                                    if (verbose){
                                      if (any(is.na(XY[, poly.id.colname]))){
                                        cat(nrow(XY[which(is.na(XY[, poly.id.colname])),]), 
                                            "points could not be mapped \n") 
                                      }
                                    } 
                                  } else {
                                    XY$sort_col <- NULL 
                                  }
                                  
                                  if (keep.data){
                                    out <- list(XY= XY, k= NA, poly.list = poly.list, 
                                                poly.id= poly.id, 
                                                poly.id.colname= poly.id.colname)
                                    return (out)
                                  } else {
                                    out <- list(XY= XY, k= NA) 
                                    return(out)
                                  }
                                  ### The original points augmented with polygon ID are returned along with the poly centers and other call information
                                  
                                }, ex = function(){
                                  data(sf.crime.2012, envir = environment())
                                  data(sf.polys, envir = environment())
                                  cat(nrow(sf.crime.2012), "rows in SF crime \n")
                                  
                                  XY.kdtree <- RapidPolygonLookup(sf.crime.2012[,c("X","Y")], poly.list= sf.polys, 
                                                                    k= 10, N= 1000, 
                                                                    poly.id= "fips", poly.id.colname= "census.block", 
                                                                    keep.data= TRUE, verbose= TRUE)
                                  
                                  XY.kdtree.DF <- XY.kdtree$XY
                                  table(XY.kdtree.DF$rank, useNA= "always")
                                  hist(XY.kdtree.DF$rank, xlab = "rank of neighbor")
                                  
                                }
  )

