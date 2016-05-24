FindPolygonInRanges <- structure(function# Use range-search to map points to polygon.
                                 ### This function searches the lat-long ranges of polygons to come up
                                 ### with a shorter list of candidates on which point.in.polygon() from the
                                 ### sp package can be applied.
                                 (poly.list, ##<< polygon list with 3 or 4 elements: data, polys, poly.centers, and possibly ranges
                                  XY, ##<< data frame containing X-Y columns
                                  poly.id = "fips", ##<< column name in poly.list$data containing the polygon identifier
                                  poly.id.colname = "census.block", ##<< desired column name in the output data frame containing the polygon identifier
                                  verbose=0 ##<< level of verbosity
                                 ){
                                   if (!("ranges"  %in% names(poly.list))) {
                                     if (verbose) print("adding missing polygon range information ...")
                                     poly.list <- AddRanges(poly.list)
                                   }
                                   #this splitting is repeated too often
                                   polys <- split(poly.list$polys[,c( "POS", "X","Y")], poly.list$polys[,c("PID")])
                                   if (!(poly.id.colname %in% colnames(XY))) XY[,poly.id.colname] <- NA
                                   jj <- which(is.na(XY[,poly.id.colname]))
                                   for (i in jj){
                                     lat <- XY[i,"Y"]
                                     lon <- XY[i,"X"]
                                     #the following should be made MUCH more efficient with binary searches:
                                     wx <- which(poly.list$ranges$x[,1]<= lon & poly.list$ranges$x[,2] >=lon)
                                     wy <- which(poly.list$ranges$y[,1]<= lat & poly.list$ranges$y[,2] >=lat)
                                     wxy <- intersect(wx,wy)
                                     if (length(wxy) == 0) {
                                       cat("point", i," is out of polygon range! \n")
                                       next
                                     }
                                     poly.candidates <- rownames(poly.list$ranges$x)[wxy]
                                     k <- 1
                                     for (pID in poly.candidates) {#need another ranking criterion here
                                       if (point.in.polygon(lon, lat, polys[[pID]][,"X"], polys[[pID]][,"Y"])>0){
                                         id <- as.character(poly.list$data[as.integer(pID), poly.id])
                                         XY[i, poly.id.colname] <- id
                                         XY[i,"rank"] <- k
                                         # if (verbose) {cat(pID, k, id, "\n")}
                                         break
                                       }
                                       k <- k+1
                                     }
                                     if (verbose) if ((i %% 5000) ==0) cat(i, " rows processed...\n")
                                   }
                                   return(XY)
                                 }, ex= function(){
                                   
                                   data(sf.crime.2012, envir = environment())
                                   data(sf.polys, envir = environment())
                                   
                                   sf.polys <- AddRanges(sf.polys)
                                   XY <- FindPolygonInRanges(sf.polys, sf.crime.2012[1:1000,], verbose=0)
                                   
                                   which(is.na(XY[,"census.block"]))
                                   table(XY$rank)

                                 }
)
                                 