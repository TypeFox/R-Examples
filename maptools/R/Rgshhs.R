# Copyright (c) 2005-2010 Roger Bivand and Karl Ove Hufthammer

Rgshhs <- function(fn, xlim=NULL, ylim=NULL, level=4, minarea=0, 
	shift=FALSE, verbose=TRUE, no.clip = FALSE, properly=FALSE,
        avoidGEOS=FALSE, checkPolygons=FALSE) {
	if (!is.character(fn)) stop("file name must be character string")
	if (length(fn) != 1L) stop("file name must be single character string")
	dolim <- FALSE
	dolim <- (!is.null(xlim) || !is.null(ylim))
	if (!is.null(xlim)) lim <- xlim
	else lim <- c(-180, 360)
	if (!is.null(ylim)) lim <- c(lim, ylim)
	else lim <- c(lim, c(-90, 90))
        storage.mode(minarea) <- "double"
	polydata <- .Call("Rgshhs", as.character(fn), as.integer(1), 
		as.logical(dolim), as.numeric(lim), as.integer(level), 
		minarea, PACKAGE="maptools")
        line <- median(polydata$line)
        if (verbose) cat("Data are", ifelse(line == 0, "polygon", "line"),
                "data\n")
	chosen_0 <- .Call("Rgshhs", as.character(fn), as.integer(2), 
		as.logical(dolim), as.numeric(lim), as.integer(level), 
		minarea, PACKAGE="maptools")
	if (dolim) clip <- .Call("Rgshhs", as.character(fn), as.integer(3), 
		as.logical(dolim), as.numeric(lim), as.integer(level), 
		minarea, PACKAGE="maptools")
	else clip <- NULL
        rgeosI <- rgeosStatus()
        if (rgeosI && !avoidGEOS) {
            # require(rgeos)
    		if (!requireNamespace("rgeos", quietly = TRUE))
				stop("rgeos spatstat required for Rgshhs")
        } else {
            stopifnot(isTRUE(gpclibPermitStatus()))
    		if (!requireNamespace("gpclib", quietly = TRUE))
				stop("gpclib spatstat required for Rgshhs")
	    	# require("gpclib")
        }
	polys <- .Call("Rgshhs", as.character(fn), as.integer(5), 
		as.logical(dolim), as.numeric(lim), as.integer(level), 
		minarea, PACKAGE="maptools")

	Antarctica <- which(polydata$area[(chosen_0+1)] > 1.3e+07 & 
		polydata$area[(chosen_0+1)] < 1.4e+07)
	if (length(Antarctica) == 1L) {
		if (verbose) cat("Polygon", which(chosen_0 == (Antarctica-1)), 
			"is Antarctica\n")
		if (verbose) cat("  area", polydata$area[Antarctica], "\n")
		if (verbose) cat("  dropping south edge to", lim[3], "\n")
		crds <- polys[[which(chosen_0 == (Antarctica-1))]]
#	    	if (verbose) print(crds[c(1,2,(nrow(crds)-5):nrow(crds)),])
		crds <- rbind(crds[1:(nrow(crds)-2),], 
			c(0, crds[(nrow(crds)-2),2]), c(0, lim[3]), 
			c(360, lim[3]), crds[1,,drop=FALSE])
#	    	if (verbose) print(crds[c(1,2,(nrow(crds)-5):nrow(crds)),])
		polys[[which(chosen_0 == (Antarctica-1))]] <- crds
	}

	if (!no.clip && dolim && any(clip == 1) && line == 0) {
	    limbb <- cbind(c(lim[1], lim[1], lim[2], lim[2], lim[1]), 
		c(lim[3], lim[4], lim[4], lim[3], lim[3]))
	    which_null <- NULL
	    opolys <- vector(mode="list", length=length(polys))
	    ic <- 1
	    if (verbose) cat("Rgshhs: clipping", sum(clip), "of", 
		length(polys), "polygons ...\n")
            if (rgeosI && !avoidGEOS) {
		limgp <- SpatialPolygons(list(Polygons(list(Polygon(limbb)),
                    ID="0")))
		for (i in seq(along=polys)) {
		    if (clip[i] == 1) {
                        tp <- SpatialPolygons(list(Polygons(list(Polygon(
                            polys[[i]])), ID="1")))
                        rp0 <- NULL
                        gI <- rgeos::gIntersection(tp, limgp)
# bug reported 120809 r-sig-geo.20.trevva
                        if (!is.null(gI) && is(gI, "SpatialCollections")) 
                            gI <- slot(gI, "polyobj")
                        if (!is.null(gI)) 
                            rp0 <- slot(gI, "polygons")[[1]]
                        rp <- NULL
                        if (!is.null(rp0)) 
                            rp <- lapply(slot(rp0, "Polygons"), slot, "coords")
			l <- length(rp)
			if (l > 0) {
		            outl <- vector(mode="list", length=l)
			    for (j in 1:l) 
				outl[[j]] <- as(rp[[j]], "matrix")
			    opolys[[ic]] <- outl
			    if (i < length(polys)) ic <- ic+1
			} else {
			    which_null <- c(which_null, i)
			    if (verbose) cat("null polygon: [[", i,
                                "]]\n", sep="");
			}
		    } else {
			opolys[[ic]] <- list(polys[[i]])
			if (i < length(polys)) ic <- ic+1
		    }
		}
		polys <- opolys[1:ic]
            } else {
		limgp <- as(limbb, "gpc.poly")
		for (i in seq(along=polys)) {
			if (clip[i] == 1) {
				tp <- as(polys[[i]], "gpc.poly")
				rp <- gpclib::intersect(tp, limgp)
				l <- length(rp@pts)
				if (l > 0) {
				    outl <- vector(mode="list", length=l)
				    for (j in 1:l) 
					outl[[j]] <- as(rp[j], "matrix")
				    opolys[[ic]] <- outl
				    if (i < length(polys)) ic <- ic+1
				} else {
					which_null <- c(which_null, i)
					if (verbose) cat("null polygon: [[",i,
					    "]]\n", sep="");
				}
			} else {
			    opolys[[ic]] <- list(polys[[i]])
			    if (i < length(polys)) ic <- ic+1
			}
		}
		polys <- opolys[1:ic]
            }
	} else {
		for (i in seq(along=polys)) polys[[i]] <- list(polys[[i]])
		which_null <- NULL
	}

	if (!is.null(which_null)) chosen_0 <- chosen_0[-which_null]
	chosen_1 <- chosen_0+1
	if (line == 0) {
	 levels <- polydata$level[chosen_1]
         if (rgeosI && !avoidGEOS) {
          ids <- polydata$id[chosen_1]
          containers <- polydata$container[chosen_1]
          ancestors <- polydata$ancestor[chosen_1]
          tl <- as.list(table(levels))
          ntl <- as.integer(names(tl))
          mntl <- match(1:4, ntl)
          l1 <- which(levels == 1L)
          if (length(l1) > 0L) {
              c1 <- which(containers == -1L)
              if (any(l1 != c1)) warning("containers and levels not coherent")
              if (!is.na(mntl[4])) {
                  wl4 <- which(levels == 4L)
                  cw4 <- containers[wl4]
                  mcw4 <- match(cw4, ids)
                  containers[wl4] <- containers[mcw4]
              }
              if (!is.na(mntl[3])) {
                  wl3 <- which(levels == 3L | levels == 4L)
                  cw3 <- containers[wl3]
                  mcw3 <- match(cw3, ids)
                  containers[wl3] <- containers[mcw3]
              }
              IDs <- ids[l1]
              if (is.na(mntl[2])) {
                  belongs <- as.list(l1)
              } else {
                  belongs <- lapply(1:length(IDs), function(i)
                      c(i, which(containers == IDs[i])))
              }
          } else {
              stop("no shoreline in selection")
          }
	  holes <- !as.logical(levels %% 2)
	  nps <- sapply(polys, length)
          n <- length(belongs)
          Srl <- vector(mode="list", length=n)
	  for (i in 1:n) {
		nParts <- length(belongs[[i]])
		srl <- NULL
		for (j in 1:nParts) {
		    this <- belongs[[i]][j]
		    for (k in 1:nps[this]) {
			crds <- polys[[this]][[k]]
			if (!identical(crds[1,], crds[nrow(crds),])) {
			    crds <- rbind(crds, crds[1,,drop=FALSE])
			    if (verbose) 
			        cat("  closing polygon", this, ":", k, "\n")
			}
			if (shift) crds[,1] <- ifelse(crds[,1] > 180, 
			    crds[,1] - 360, crds[,1])
			if (checkPolygons) {
                            jres <- list(Polygon(crds))
                        } else {
                            jres <- list(Polygon(crds, hole=holes[this]))
                        }
			srl <- c(srl, jres)
		    }
		}
                pls0 <- Polygons(srl, ID=IDs[i])
		if (checkPolygons) {
                    Srl[[i]] <- checkPolygonsGEOS(pls0, properly=properly)
                } else {
                    Srl[[i]] <- pls0
                }
	  }
	  res <- as.SpatialPolygons.PolygonsList(Srl, 
		proj4string=CRS("+proj=longlat +datum=WGS84"))
          polydata <- data.frame(polydata)[chosen_1,]

	  return(list(polydata=polydata, belongs=belongs, SP=res))
         } else {
	  belongs <- matrix(1:length(chosen_1), ncol=1)
#	  belonged_to <- as.numeric(rep(NA, length(chosen_1)))


	  if (level > 1 && any(levels > 1)) {
	    if (verbose) {
		cat("Rgshhs: assigning enclosed polygons to their enclosers\n")
		cat("  level tallies:\n")
		print(table(levels))
		cat("...\n")
	    }
	    mlevel <- as.integer(max(levels))
	    belongs <- matrix(rep(1:length(chosen_1), mlevel), ncol=mlevel)
	    first_time <- TRUE
	    for (il in mlevel:2) {
		w_il <- which(levels == il)
		w_il_1 <- which(levels == (il-1))
		if (length(w_il) > 0L) {
			if (length(w_il_1) == 1L) {
			    belongs[w_il, (il-1)] <- w_il_1
			    if (!first_time) {
				prom <- which(!is.na(match(belongs[,il], w_il)))
				belongs[prom, (il-1)] <- rep(w_il_1, 
				    length(prom))
			    }
			    first_time <- FALSE
			} else {
			    l_1 <- vector(mode="list", length=length(w_il_1))
			    for (i in 1:length(w_il_1)) {
				ii <- w_il_1[i]
				lp1 <- as(polys[[ii]][[1]], "gpc.poly")
				if (length(polys[[ii]]) > 1L) {
				    for (j in 2:length(polys[[ii]])) {
					lpj <- as(polys[[ii]][[j]], "gpc.poly")
					lp1 <- gpclib::append.poly(lp1, lpj)
				    }
				}
				l_1[[i]] <- lp1
			    }
			    for (i in 1:length(w_il)) {
				ii <- w_il[i]
				lp1 <- as(polys[[ii]][[1]], "gpc.poly")
				if (length(polys[[ii]]) > 1L) {
				    for (j in 2:length(polys[[ii]])) {
					lpj <- as(polys[[ii]][[j]], "gpc.poly")
					lp1 <- gpclib::append.poly(lp1, lpj)
				    }
				}
				for (j in 1:length(l_1)) {
				    tp <- gpclib::intersect(l_1[[j]], lp1)
				    if (length(tp@pts) > 0L) {
					belongs[w_il[i], (il-1)] <- w_il_1[j]
			    		if (!first_time) {
					    prom <- which(!is.na(match(
						belongs[,il], w_il[i])))
					    belongs[prom, (il-1)] <- w_il_1[j]
					}
					break
				    }
				}
			    }
			    first_time <- FALSE
			}
		}
	    }
	  }

	  if (verbose) cat("Rgshhs: constructing SpatialPolygons ...\n")
	  holes <- !as.logical(levels %% 2)
	  nps <- sapply(polys, length)
	  IDs <- polydata$id[chosen_1[belongs[,1]]]
	  tab <- table(factor(IDs))
	  n <- length(tab)
	  IDss <- names(tab)
	  reg <- match(IDs, IDss)
	  new_belongs <- lapply(1:n, function(x) which(x == reg))
	  Srl <- vector(mode="list", length=n)
	  for (i in 1:n) {
		nParts <- length(new_belongs[[i]])
		srl <- NULL
		for (j in 1:nParts) {
		    this <- new_belongs[[i]][j]
		    for (k in 1:nps[this]) {
			crds <- polys[[this]][[k]]
			if (!identical(crds[1,], crds[nrow(crds),])) {
			    crds <- rbind(crds, crds[1,,drop=FALSE])
			    if (verbose) 
			        cat("  closing polygon", this, ":", k, "\n")
			}
			if (shift) crds[,1] <- ifelse(crds[,1] > 180, 
			    crds[,1] - 360, crds[,1])
			jres <- list(Polygon(crds, hole=holes[this]))
			srl <- c(srl, jres)
		    }
		}
		Srl[[i]] <- Polygons(srl, ID=IDss[i])
	  }
	  res <- as.SpatialPolygons.PolygonsList(Srl, 
		proj4string=CRS("+proj=longlat +datum=WGS84"))
	  list(polydata=data.frame(polydata)[chosen_1,], belongs=belongs,
		new_belongs=new_belongs, SP=res)
         }
	} else {
	  Sll <- lapply(1:length(polys), function(i) {
              ID <- as.character(i)
              crds <- polys[[i]][[1]]
		if (shift) crds[,1] <- ifelse(crds[,1] > 180, 
		    crds[,1] - 360, crds[,1])
              Ln <- Line(crds)
              Lines(list(Ln), ID=ID)
            })
          res <- SpatialLines(Sll, 
            proj4string=CRS("+proj=longlat +datum=WGS84"))
	  list(SP=res)
	}
}

# contributed 101018 by Karl Ove Hufthammer

getRgshhsMap = function (fn = system.file("share/gshhs_c.b",
 package = "maptools"), xlim, ylim, level = 1, shift = TRUE,
 verbose = TRUE, no.clip = FALSE, properly=FALSE, avoidGEOS=FALSE, checkPolygons=FALSE) 
{
    # First try fetching the map directly, possibly with negative coordinates.
    # Note that some polygons with longitude < 0 use negative coordinates 
    # (e.g., Great Britain), and some use positve coordinates (e.g., Ireland).
    #    
    # (Must use 'try' here, because for example xlim=c(-40,-10)
    # results in an error, while xlim=c(-40,-5) does not.)
    map1 = try(Rgshhs(fn, xlim = xlim, ylim = ylim, shift = shift, 
                    level = level, verbose=verbose, no.clip = no.clip,
                    properly=properly, avoidGEOS=avoidGEOS,
                    checkPolygons=checkPolygons)$SP)
    
    # Now try fetching the same area using positive coordinates.
    xl.west = (xlim + 360)%%360
    if (xl.west[2] < xl.west[1])
        xl.west[2] = 360
    map2 = Rgshhs(fn, xlim = xl.west, ylim = ylim, shift = shift, 
            level = level, verbose=verbose, no.clip = no.clip,
            properly=properly, avoidGEOS=avoidGEOS,
            checkPolygons=checkPolygons)$SP
    
    # If there where no polygons with negative coordinates, just
    # use the positive coordinates.
    if (class(map1) == "try-error") 
        map.union = map2 else { # Else merge the two maps into one.
        
        # First store the original polygon IDs in data frames.
        df1 = data.frame(polyID = row.names(map1), stringsAsFactors=FALSE)
        row.names(df1) = df1$polyID
        map1.spdf = SpatialPolygonsDataFrame(map1, df1)
        
        df2 = data.frame(polyID = row.names(map2), stringsAsFactors=FALSE)
        row.names(df2) = df2$polyID
        map2.spdf = SpatialPolygonsDataFrame(map2, df2)
        
        # Generate new polygon IDs to avoid duplicate IDs when
        # rbinding the two maps.
        row.names(map1.spdf) = as.character(seq_along(map1@polygons))
        row.names(map2.spdf) = as.character(length(map1@polygons) + 
                        seq_along(map2@polygons))
        map.merged = rbind(map1.spdf, map2.spdf)
        
        # Finally, combine all the polygons, using the
        # original polyon IDs.
        map.union = unionSpatialPolygons(map.merged, map.merged$polyID)
    }
    map.union
}


