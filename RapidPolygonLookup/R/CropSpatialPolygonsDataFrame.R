CropSpatialPolygonsDataFrame <-structure(function# Crop polygons to bounding box and adds polygon centers
                                         ### This function serves three purposes:
                                         ### (i)   changes the (complicated) data structure of a spatial polygon (from the sp package) to a format which is aligned with the (simpler) PBSmapping polygon format.
                                         ### (ii)  clips/crops the polygons to a pre specified bounding box
                                         ### (iii) computes and adds the polygon centers for each polygon 
                                         (
                                           x,##<< object of class SpatialPolygonsDataFrame
                                           bb= NULL,##<< bounding box to crop the polygons
                                           verbose= 0 ##<< level of verbosity
                                         ){
                                           N <- nrow(x@data)
                                           
                                           data <- list()
                                           polys.to.keep <- list()
                                           poly.centers <- list()
                                           
                                           if (is.null(bb)){
                                             bb <- data.frame(t(bbox(x)))
                                             names(bb) <- c("X", "Y")
                                           }
                                           
                                           k=1
                                           if (verbose >1) browser()
                                           for (i in 1:N) {
                                             polys = x@polygons[[i]]@Polygons[[1]]@coords
                                             M = nrow(polys)
                                             polys <- data.frame(PID=rep(k, M), POS= 1:M, X= polys[,1], Y= polys[,2])
                                             new.polys <- clipPolys(polys, xlim= bb$X, ylim= bb$Y)
                                             if (!is.null(new.polys)) {
                                               polys.to.keep[[k]] <- new.polys
                                               data[[k]] <- x@data[i,,drop=F]
                                               poly.centers[[k]] = calcCentroid(new.polys, rollup = 3)
                                               k <- k+1
                                             }
                                           }
                                           data = do.call("rbind", data)
                                           polys.to.keep = do.call("rbind", polys.to.keep)
                                           poly.centers = do.call("rbind", poly.centers)
                                           
                                           result <- list(data = data, polys = polys.to.keep, poly.centers = poly.centers)
                                           return(result)
                                           ### New list with seperate entries for data, polys, and poly centers
                                         }, ex = function(){                                           
                                           # San Francisco:
                                           data(california.tract10, envir = environment())
                                           sf.polys <- CropSpatialPolygonsDataFrame(x= california.tract10, 
                                                                                  bb= data.frame(X=c(-122.5132, -122.37), 
                                                                                                 Y= c(37.70760, 37.81849)))
                                           
                                         })

