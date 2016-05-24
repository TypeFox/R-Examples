SearchForPolygon <- structure(function# Use kd-trees to search the nearest neighbour polygons for a given set of points
                              ### This function uses the nn2() function from the RANN package to come up
                              ### with a shorter list of candidates on which point.in.polygon() from the 
                              ### sp package can be applied.
                              (poly.list, ##<< polygon list with 3-4 elements: poly.centers, data, polys and possibly ranges
                               XY, ##<< data frame containing X-Y columns to assign polygons to
                               k, ##<< maximum number of nearest neighbours to compute. The default value is set to 10.
                               poly.id, ##<< column name in poly.list$data containing the polygon identifier
                               poly.id.colname, ##<< desired column name in the output data frame containing the polygon identifier
                               verbose= 0 ##<< level of verbosity
                               ){
                                nearest <- nn2(poly.list$poly.centers[,c("X","Y")] , query= XY[,c("X","Y")] , k= k)
                                K <- ncol(nearest$nn.idx)#should be equal to k !
                                if (verbose>1) browser()
                                polys <- split(poly.list$polys[,c( "POS", "X","Y")], poly.list$polys[,c("PID")])
                                N <- nrow(XY)
                                XY$sort_col <- 1:nrow(XY)
                                XY <- lapply(1:N, 
                                             function(i){
                                               for (j in 1:K){
                                                 try({  
                                                   pID <- nearest$nn.idx[i,j]
                                                   if (point.in.polygon(XY[i, "X"], XY[i, "Y"], polys[[pID]][,"X"], polys[[pID]][,"Y"])>0){
                                                     id <- as.character(poly.list$data[pID, poly.id])
                                                     rank <- j
                                                     break
                                                   } else{
                                                     id <- NA
                                                     rank <- NA
                                                   } 
                                                 })
                                               }
                                               if (verbose) if ((i %% 5000) ==0) cat(i, " rows processed...\n")
                                               out <- data.frame(X= XY[i, "X"], Y= XY[i, "Y"], id= id, rank= rank, 
                                                                 sort_col = XY[i, "sort_col"], stringsAsFactors= FALSE)
                                               names(out)[3] <- poly.id.colname
                                               return(out)
                                             })
                                tmp <- names(XY[[1]])
                                XY <- data.frame(matrix(unlist(XY), nrow= length(XY), byrow= TRUE), 
                                                 stringsAsFactors= FALSE)
                                names(XY) <- tmp
                                XY$X <- as.numeric(XY$X)
                                XY$Y <- as.numeric(XY$Y)
                                XY$rank <- as.numeric(XY$rank)
                                XY$sort_col <- as.numeric(XY$sort_col)
                                XY
                                ### Returns data frame with identified polygon and nearest neighbour rank
                              }, ex= function(){
                                data(sf.crime.2012, envir = environment())
                                data(sf.polys, envir = environment())
                                XY.polys <- SearchForPolygon(poly.list= sf.polys, XY= sf.crime.2012[1:1000, ], k= 10,
                                                             poly.id= "fips", poly.id.colname= "census.block",
                                                             verbose= TRUE)
                                
                              }
  )
  
  
  