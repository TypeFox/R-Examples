.adjFormat <- function(x) {
  x <- x[, 1:3]
  names(x) <- c("identifier", "XCOOR", "YCOOR")
  return(x)
}
.BarChartPoly <- function(x, plotout = F, verbose = FALSE, ...) {
    if (!class(x) == "spgeoOUT" && !class(x) == "spgeoH") {
        stop("this function is only defined for class \"spgeoOUT\"")
    }
    if (plotout == TRUE) {
        par(ask = F)
    }
    if (plotout == FALSE) {
      par(ask = T)
    }
    liste <- names(x$spec_table)
    leng <- length(liste)
    if (length(names(x$spec_table)) == 0) {
        cat("No point fell in any polygon")
    } else {
        for (i in 2:leng) {
            subs <- subset(x$spec_table, x$spec_table[, i] > 0)
            datsubs <- subs[order(subs[, i]), ]
            if (dim(subs)[1] == 0) {
                plot(1:10, 1:10, type = "n", xlab = "", ylab = "Number of occurences")
                text(3, 6, labels = "No species occurred in this polygon.", adj = 0)
                title(liste[i])
            } else {
                barplot(datsubs[, i], names.arg = datsubs$identifier, 
                        las = 2, ylab = "Number of occurences", cex.names = 0.7,
                        ylim = c(0, (max(datsubs[, i]) + max(datsubs[, i])/10)))  #, ...)
                box("plot")
                title(liste[i])
            }
        }
    }
    par(ask = F)
    
}
.BarChartSpec <- function(x, mode = c("percent", "total"), plotout = FALSE, verbose = FALSE, ...) {
    match.arg(mode)
    if (!class(x) == "spgeoOUT" && !class(x) == "spgeoH") {
        stop("this function is only defined for class \"spgeoOUT\"")
    }
    if (length(x$spec_table) == 0) {
        cat("No point was found inside the given polygons")
    } else {
        if (plotout == FALSE) {
            par(ask = T)
        }
        if (mode[1] == "total") {
            liste <- x$spec_table$identifier
            leng <- length(liste)
            #par(mar = c(10, 4, 3, 3))
            for (i in 1:leng) {
                if (verbose == TRUE) {
                  cat(paste("Creating barchart for species ", i, "/", leng, ": ", liste[i], "\n", sep = ""))
                }
                spsub <- as.matrix(subset(x$spec_table, x$spec_table$identifier == liste[i])[, 2:dim(x$spec_table)[2]])
                if (sum(spsub) > 0) {
                  barplot(spsub, las = 2, ylim = c(0, (max(spsub) + max(spsub)/10)), ylab = "Number of occurrences", ...)
                  title(liste[i])
                  box("plot")
                }
            }
        }
        if (mode[1] == "percent") {
            percent <- x$spec_table[, -1]
            anzpoly <- length(names(x$spec_table)[-1])
            if (anzpoly > 1) {
                percent2 <- percent/rowSums(percent) * 100
            } else {
                percent2 <- percent/sum(percent) * 100
            }
            percent2[percent2 == "NaN"] <- 0
            percent2 <- data.frame(identifier = x$spec_table[, 1], percent2)
            
            liste <- x$spec_table$identifier
            leng <- length(liste)
            leng2 <- length(colnames(percent2))
            #par(mar = c(10, 4, 3, 3))
            for (i in 1:leng) {
                if (verbose == TRUE) {
                  cat(paste("Creating barchart for species ", i, "/", leng, ": ", liste[i], "\n", sep = ""))
                }
                if (anzpoly > 1) {
                  spsub <- as.matrix(subset(percent2, percent2$identifier == liste[i])[, 2:leng2])
                } else {
                  spsub <- as.matrix(percent2[percent2$identifier == liste[i], ][, 2:leng2])
                  names(spsub) <- names(x$spec_table)[-1]
                }
                if (sum(spsub) > 0) {
                  barplot(spsub, las = 2, ylim = c(0, (max(spsub) + max(spsub)/10)), ylab = "Percent of occurrences", names.arg = names(spsub), 
                    ...)
                  title(liste[i])
                  box("plot")
                }
            }
        }
        par(ask = F)
    }
}
.CoExClassH <- function(x, verbose = FALSE) {
    dat <- x
    if (length(dim(dat)) == 0) {
        coemat <- "NULL"
    } else {
        if (!is.data.frame(x)) {
            stop("function only defined for class \"data.frame\"")
        }
        if ("identifier" %in% names(dat) == F) {
            if (T %in% sapply(dat, is.factor)) {
                id <- sapply(dat, is.factor)
                old <- names(dat)[id == T]
                names(dat)[id == T] <- "identifier"
                warning(paste("no species identifier found in input object. \n", "Column <", old, "> was used as identifier", 
                  sep = ""))
            }
            if (T %in% sapply(dat, is.character)) {
                id <- sapply(dat, character)
                old <- names(dat)[id == T]
                names(dat)[id == T] <- "identifier"
                warning(paste("no species identifier found in input object. \n", "Column <", old, "> was used as identifier", 
                  sep = ""))
            }
        }
        spnum <- length(dat$identifier)
        numpol <- length(names(dat))
        coemat <- data.frame(matrix(NA, nrow = spnum, ncol = spnum))
        for (j in 1:spnum) {
            if (verbose == TRUE) {
                cat(paste("Calculate coexistence pattern for species: ", j, "/", spnum, " ", dat$identifier[j], "\n", sep = ""))
            }
            sco <- data.frame(dat$identifier)
            for (i in 2:length(names(dat))) {
                if (dat[j, i] == 0) {
                  poly <- rep(0, spnum)
                  sco <- cbind(sco, poly)
                }
                if (dat[j, i] > 0) {
                  scoh <- dat[, i]
                  if (numpol > 2) {
                    totocc <- rowSums(dat[j, -1])
                  } else {
                    totocc <- dat[j, -1]
                  }
                  for (k in 1:length(scoh)) if (scoh[k] > 0) {
                    scoh[k] <- dat[j, i]/totocc * 100
                  } else {
                    scoh[k] <- 0
                  }
                  sco <- cbind(sco, scoh)
                }
            }
            if (numpol > 2) {
                coex <- rowSums(sco[, -1])
                coemat[j, ] <- coex
            } else {
                coex <- sco[, -1]
                coemat[j, ] <- coex
            }
        }
        coemat <- cbind(dat$identifier, coemat)
        names(coemat) <- c("identifier", as.character(dat$identifier))
    }
    return(coemat)
}
.ConvHull <- function(x){
  conv.hull <- chull(x$XCOOR, x$YCOOR)
  dat2 <- x[conv.hull, ]
  dat2 <- rbind(dat2[, c(2, 3)], dat2[1, c(2, 3)])
  poly <- SpatialPolygons(list(Polygons(list(Polygon(dat2)), ID = paste(x[1, 1], "_convhull", sep = ""))), proj4string = CRS("+proj=longlat +datum=WGS84"))
  return(poly)
}

.ConvertPoly <- function(x) {
  x <- read.table(x, sep = "\t")
  
  out2 <- vector()
  
  for (j in 1:dim(x)[1]) {
    aa <- as.character(x[j, ])
    ff <- t(aa)
    bb <- unlist(strsplit(ff[1], split = ":"))
    bb <- c(bb[1], unlist(strsplit(bb[2], split = " ")))
    
    out <- c(1, 1, 1)
    
    for (i in 2:length(bb)) {
      dd <- c(bb[1], unlist(strsplit(as.character(bb[i]), split = ",")))
      out <- rbind(out, dd)
    }
    out2 <- rbind(out2, out[-c(1, 2), ])
  }
  
  colnames(out2) <- c("identifier", "XCOOR", "YCOOR")
  rownames(out2) <- NULL
  out2 <- as.data.frame(out2)
  return(out2)
}

.Cord2Polygon <- function(x) {
    if (is.character(x)) {
        tt <- read.table(x, sep = "\t")
        if (dim(tt)[2] != 3) {
            stop(paste("wrong input format: \n", "Inputobject must be a tab-delimited text file or a data.frame with three columns", 
                sep = ""))
        }
        if (!is.numeric(tt[, 2]) || !is.numeric(tt[, 3])) {
            stop(paste("wrong input format: \n", "Input coordinates (columns 2 and 3) must be numeric", sep = ""))
        }
        if (!is.character(tt[, 1]) && !is.factor(tt[, 1])) {
            warning("input identifier (column 1) should be a string or a factor")
        }
        names(tt) <- c("identifier", "lon", "lat")
        liste <- levels(tt$identifier)
        col <- list()
        for (i in 1:length(liste)) {
            pp <- subset(tt, tt$identifier == liste[i])[, c(2, 3)]
            pp <- Polygon(pp)
            po <- Polygons(list(pp), ID = liste[i])
            col[[i]] <- po
        }
        polys <- SpatialPolygons(col, proj4string = CRS("+proj=longlat +datum=WGS84"))
    } else {
        tt <- x
        if (dim(tt)[2] != 3) {
            stop(paste("wrong input format: \n", "Inputobject must be a tab-delimited text file or a data.frame with three columns", 
                sep = ""))
        }
        if (!is.numeric(tt[, 2]) || !is.numeric(tt[, 3])) {
            stop(paste("wrong input format: \n", "Input coordinates (columns 2 and 3) must be numeric", sep = ""))
        }
        if (!is.character(tt[, 1]) && !is.factor(tt[, 1])) {
            warning("input identifier (column 1) should be a string or a factor")
        }
        names(tt) <- c("identifier", "lon", "lat")
        liste <- levels(tt$identifier)
        col <- list()
        for (i in 1:length(liste)) {
            pp <- subset(tt, tt$identifier == liste[i])[, c(2, 3)]
            pp <- Polygon(pp)
            po <- Polygons(list(pp), ID = liste[i])
            col[[i]] <- po
        }
        polys <- SpatialPolygons(col, proj4string = CRS("+proj=longlat +datum=WGS84"))
    }
    return(polys)
}
.eoo <- function(x) {
  if (!requireNamespace("geosphere", quietly = TRUE)) {
    stop("geosphere needed for this function to work. Please install it.",
         call. = FALSE)
  }  
    conv.hull <- chull(x$XCOOR, x$YCOOR)
    dat2 <- x[conv.hull, ]
    dat2 <- rbind(dat2[, c(2, 3)], dat2[1, c(2, 3)])
    poly <- SpatialPolygons(list(Polygons(list(Polygon(dat2)), ID = paste(x[1, 1], "_convhull", sep = ""))), proj4string = CRS("+proj=longlat +datum=WGS84"))
    area <- round(geosphere::areaPolygon(poly)/(1000 * 1000), 0)
    return(area)
}
.getEle <- function(x) {
  ele <- try(getData("SRTM", lon = round(as.numeric(x[2]), 2), lat = round(as.numeric(x[3]), 2)))
  if (class(ele) == "try-error") {
    elevation <- "NA"
  } else {
    if (!is.na(extract(ele, SpatialPoints(data.frame(round(as.numeric(x[2]), 2), round(as.numeric(x[3]), 2)))))) {
      elevation <- extract(ele, data.frame(round(as.numeric(x[2]), 2), round(as.numeric(x[3]), 2)))
    } else {
      elevation <- "NA"
    }
  }
  return(elevation)
}
.GetPythonIn <- function(inpt) {
    
    coord <- read.table(inpt[1], header = T, sep = "\t")
    idi <- coord[, 1]
    coords <- coord[, c(2, 3)]
    
    polyg <- read.table(inpt[2], header = T, sep = "\t")
    poly <- .Cord2Polygon(polyg)
    
    samtab <- read.table(inpt[3], header = T, sep = "\t")
    
    spectab <- read.table(inpt[4], header = T, sep = "\t")
    names(spectab)[1] <- "identifier"
    
    polytab <- .SpPerPolH(spectab)
    
    nc <- subset(samtab, is.na(samtab$homepolygon))
    identifier <- idi[as.numeric(rownames(nc))]
    bb <- coords[as.numeric(rownames(nc)), ]
    noclass <- data.frame(identifier, bb)
    
    
    outo <- list(identifier_in = idi, species_coordinates_in = coords, polygons = poly, sample_table = samtab, spec_table = spectab, 
        polygon_table = polytab, not_classified_samples = noclass, coexistence_classified = "NA")
    class(outo) <- "spgeoOUT"
    return(outo)
}
.HeatPlotCoEx <- function(x, verbose = FALSE, ...) {
  if (class(x) == "spgeoOUT") {
    dat <- x$coexistence_classified
  } else {
    dat <- x
  }
  if (dim(dat)[1] > 40) {
    warning("more than 40 species in coexistence matrix. Plot might be unreadable")
  }
  if (class(dat) != "data.frame") {
    stop("wrong input format. Input must be a \"data.frame\"")
  }
  if (dim(dat)[2] != (dim(dat)[1] + 1)) {
    warning("suspicicous data dimensions, check input file")
  }
  ymax <- dim(dat)[1]
  xmax <- dim(dat)[2]
  colo <- rev(heat.colors(10))
  numer <- rev(1:ymax)
  
  layout(matrix(c(rep(1, 9), 2), ncol = 1, nrow = 10))
  par(mar = c(0, 10, 10, 0))
  plot(0, xlim = c(0, xmax - 1), ylim = c(0, ymax), type = "n", axes = F, xlab = "", ylab = "")
  for (j in 2:xmax) {
    for (i in 1:ymax) {
      if (i == (j - 1)) {
        rect(j - 2, numer[i] - 1, j - 1, numer[i], col = "black")
      } else {
        ind <- round(dat[i, j]/10, 0)
        if (ind == 0) {
          rect(j - 2, numer[i] - 1, j - 1, numer[i], col = "white")
        } else {
          rect(j - 2, numer[i] - 1, j - 1, numer[i], col = colo[ind])
        }
      }
    }
  }
  axis(side = 3, at = seq(0.5, (xmax - 1.5)), labels = colnames(dat)[-1], las = 2, cex.axis = 0.7, pos = ymax)
  axis(2, at = seq(0.5, ymax), labels = rev(dat$identifier), las = 2, cex.axis = 0.7, pos = 0)
  title("Species co-occurrence", line = 9)
  
  par(mar = c(0.5, 10, 0, 0))
  plot(c(1, 59), c(1, 12), type = "n", axes = F, ylab = "", xlab = "")
  text(c(13, 13), c(10, 7), c("0%", "10%"))
  text(c(20, 20), c(10, 7), c("20%", "30%"))
  text(c(27, 27), c(10, 7), c("40%", "50%"))
  text(c(34, 34), c(10, 7), c("60%", "70%"))
  text(c(41, 41), c(10, 7), c("80%", "90%"))
  text(c(48), 10, "100%")
  rect(c(9, 9, 16, 16, 23, 23, 30, 30, 37, 37, 44), c(rep(c(10.7, 7.7), 5), 10.7), 
       c(11, 11, 18, 18, 25, 25, 32, 32, 39, 39, 46), c(rep(c(8.7, 
                                                                                                                                          5.7), 5), 8.7), col = c("white", colo))
  rect(7, 5, 51, 12)
  layout(matrix(1, 1, 1))
} 
.MapAll <- function(x, polyg, moreborders = FALSE, verbose = FALSE, ...) {
    # data('wrld_simpl', envir = environment())
    if (class(x) == "spgeoOUT") {
        xmax <- min(max(x$species_coordinates_in[, 1]) + 2, 180)
        xmin <- max(min(x$species_coordinates_in[, 1]) - 2, -180)
        ymax <- min(max(x$species_coordinates_in[, 2]) + 2, 90)
        ymin <- max(min(x$species_coordinates_in[, 2]) - 2, -90)
        difx <- sqrt(xmax^2 + xmin^2)
        dify <- sqrt(ymax^2 + ymin^2)
        if (difx > 90) {
            xmax <- min(xmax + 10, 180)
            xmin <- max(xmin - 10, -180)
            ymax <- min(ymax + 10, 90)
            ymin <- max(ymin - 10, -90)
        }
        if (verbose == TRUE) {
            warning("creating map of all samples")
        }
        map("world", xlim = c(xmin, xmax), ylim = c(ymin, ymax))
        axis(1)
        axis(2)
        box("plot")
        title("All samples")
        # if (moreborders == T) {plot(wrld_simpl, add = T)}
        if (verbose == TRUE) {
            warning("adding polygons")
        }
        plot(x$polygons, col = "grey60", border = "grey40", add = T, ...)
        if (verbose == TRUE) {
            warning("adding sample points")
        }
        points(x$species_coordinates_in[, 1], x$species_coordinates_in[, 2], cex = 0.7, pch = 3, col = "blue", ...)
    }
    if (class(x) == "matrix" || class(x) == "data.frame") {
        if (!is.numeric(x[, 1]) || !is.numeric(x[, 2])) {
            stop(paste("wrong input format:\n", "Point input must be a \"matrix\" or \"data.frame\" with 2 columns.\n", "Column order must be lon - lat", 
                sep = ""))
        }
        if (class(polyg) != "SpatialPolygons") {
            warning("to plot polygons, polyg must be of class \"SpatialPolygons\"")
        }
        x <- as.data.frame(x)
        nums <- sapply(x, is.numeric)
        x <- x[, nums]
        xmax <- min(max(x[, 2]) + 2, 180)
        xmin <- max(min(x[, 2]) - 2, -180)
        ymax <- min(max(x[, 1]) + 2, 90)
        ymin <- max(min(x[, 1]) - 2, -90)
        if (ymax > 92 || ymin < -92) {
            warning("column order must be lon-lat, not lat - lon. Please check")
        }
        map("world", xlim = c(xmin, xmax), ylim = c(ymin, ymax))
        axis(1)
        axis(2)
        title("All samples")
        box("plot")
        # if (moreborders == T) {plot(wrld_simpl, add = T, ...)}
        if (class(polyg == "list")) 
            
        plot(polyg, col = "grey60", add = T, ...)
        
        points(x[, 2], x[, 1], cex = 0.5, pch = 3, col = "blue", ...)
        dat <- data.frame(x$not_classified_samples)
        points(dat$XCOOR, dat$YCOOR, cex = 0.5, pch = 3, col = "red", ...)
    }
}
.MapPerPoly <- function(x, areanames = NULL, plotout = FALSE) {
    if (!class(x) == "spgeoOUT") {
        stop("this function is only defined for class \"spgeoOUT\"")
    }
    dum <- x$polygons
    
    if (class(dum) == "SpatialPolygonsDataFrame") {
      if (length(areanames) == 0) {
        areanames <- x$areanam
      }
        len <- length(unique(x$polygons@data[, areanames]))
    } else {
        len <- length(names(dum))
    }
    for (i in 1:len) {
        if (class(dum) == "SpatialPolygonsDataFrame") {
            chopo <- unique(x$polygons@data[, areanames])[i]
            
            xmax <- min(max(bbox(subset(x$polygons, x$polygons@data[, areanames] == unique(x$polygons@data[, areanames])[i]))[1, 
                2]) + 5, 180)
            xmin <- max(min(bbox(subset(x$polygons, x$polygons@data[, areanames] == unique(x$polygons@data[, areanames])[i]))[1, 
                1]) - 5, -180)
            ymax <- min(max(bbox(subset(x$polygons, x$polygons@data[, areanames] == unique(x$polygons@data[, areanames])[i]))[2, 
                2]) + 5, 90)
            ymin <- max(min(bbox(subset(x$polygons, x$polygons@data[, areanames] == unique(x$polygons@data[, areanames])[i]))[2, 
                1]) - 5, -90)
        } else {
            chopo <- names(dum)[i]
            
            xmax <- min(max(bbox(x$polygons[i])[1, 2]) + 5, 180)
            xmin <- max(min(bbox(x$polygons[i])[1, 1]) - 5, -180)
            ymax <- min(max(bbox(x$polygons[i])[2, 2]) + 5, 90)
            ymin <- max(min(bbox(x$polygons[i])[2, 1]) - 5, -90)
        }
        po <- data.frame(x$sample_table, x$species_coordinates_in)
        subpo <- subset(po, as.character(po$homepolygon) == as.character(chopo))
        
        subpo <- subpo[order(subpo$identifier), ]
        
        liste <- unique(subpo$identifier)
        leng <- length(liste)
        
        rain <- rainbow(leng)
        ypos <- vector(length = leng)
        yled <- (ymax - ymin) * 0.025
        for (k in 1:leng) {
            ypos[k] <- ymax - yled * k
        }
        
        layout(matrix(c(1, 1, 1, 1, 1, 2, 2), ncol = 7, nrow = 1))
        par(mar = c(3, 3, 3, 0))
        te <- try(map("world", xlim = c(xmin, xmax), ylim = c(ymin, ymax)), silent = T)
        if (class(te) == "try-error") {
            map("world")
        }
        axis(1)
        axis(2)
        box("plot")
        title(chopo)
        if (class(dum) == "SpatialPolygonsDataFrame") {
            plot(subset(x$polygons, x$polygons@data[, areanames] == unique(x$polygons@data[, areanames])[i]), col = "grey60", 
                add = T)
        } else {
            plot(x$polygons[i], col = "grey60", add = T)
        }
        for (j in 1:leng) {
            subsub <- subset(subpo, subpo$identifier == liste[j])
            points(subsub[, 3], subsub[, 4], cex = 1, pch = 3, col = rain[j])
        }
        par(mar = c(3, 0, 3, 0), ask = F)
        plot(c(1, 50), c(1, 50), type = "n", axes = F)
        if (leng == 0) {
            yset <- 25
            xset <- 1
        }
        if (leng == 1) {
            yset <- 25
            xset <- rep(4, leng)
        }
        if (leng > 1) {
            yset <- rev(sort(c(seq(25, 25 + max(ceiling(leng/2) - 1, 0)), seq(24, 24 - leng/2 + 1))))
            xset <- rep(4, leng)
        }
        points(xset - 2, yset, pch = 3, col = rain)
        if (leng == 0) {
            text(xset, yset, labels = "No species found in this polygon", adj = 0)
        } else {
            text(xset, yset, labels = liste, adj = 0, xpd = T)
            rect(min(xset) - 4, min(yset) - 1, 50 + 1, max(yset) + 1, xpd = T)
        }
        
        if (plotout == FALSE) {
            par(ask = T)
        }
    }
    par(ask = F)
    layout(matrix(1,1,1))
}
.MapPerSpecies <- function(x, moreborders = FALSE, plotout = FALSE, verbose = FALSE, ...) {
    if (!class(x) == "spgeoOUT") {
        stop("this function is only defined for class \"spgeoOUT\"")
    }
    # if (moreborders == T) {data('wrld_simpl', envir = environment())}
    layout(matrix(1, ncol = 1, nrow = 1))
    if (plotout == FALSE) {
        par(ask = T)
    }
    dat <- data.frame(x$sample_table, x$species_coordinates_in)
    names(dat) <- c("identifier", "homepolygon", "XCOOR", "YCOOR")
    liste <- unique(dat$identifier)
    
    
    for (i in 1:length(liste)) {
#         if (verbose == TRUE) {
#             cat(paste("Mapping species: ", i, "/", length(liste), ": ", liste[i], "\n", sep = ""))
#         }
        kk <- subset(dat, dat$identifier == liste[i])
        
        inside <- kk[!is.na(kk$homepolygon), ]
        outside <- kk[is.na(kk$homepolygon), ]
        
        xmax <- min(max(dat$XCOOR) + 2, 180)
        xmin <- max(min(dat$XCOOR) - 2, -180)
        ymax <- min(max(dat$YCOOR) + 2, 90)
        ymin <- max(min(dat$YCOOR) - 2, -90)
        
        map("world", xlim = c(xmin, xmax), ylim = c(ymin, ymax))
        axis(1)
        axis(2)
        title(liste[i])
        # if (moreborders == T) {plot(wrld_simpl, add = T)}
        plot(x$polygons, col = "grey60", add = T)
        
        if (length(inside) > 0) {
            points(inside$XCOOR, inside$YCOOR, cex = 0.7, pch = 3, col = "blue")
        }
        if (length(outside) > 0) {
            points(outside$XCOOR, outside$YCOOR, cex = 0.7, pch = 3, col = "red")
        }
        box("plot")
    }
    par(ask = F)
    layout(matrix(1,1,1))
}
.MapUnclassified <- function(x, moreborders = FALSE, verbose = FALSE, ...) {
    if (!class(x) == "spgeoOUT") {
        stop("This function is only defined for class \"spgeoOUT\"")
    }
    dat <- data.frame(x$not_classified_samples)
    # if (moreborders == T) {data('wrld_simpl', envir = environment())}
    if (dim(dat)[1] == 0) {
        plot(c(1:20), c(1:20), type = "n", axes = F, xlab = "", ylab = "")
        text(10, 10, labels = paste("All points fell into the polygons and were classified.\n", "No unclassified points", sep = ""))
    } else {
        xmax <- min(max(dat$XCOOR) + 2, 180)
        xmin <- max(min(dat$XCOOR) - 2, -180)
        ymax <- min(max(dat$YCOOR) + 2, 90)
        ymin <- max(min(dat$YCOOR) - 2, -90)
        
        map("world", xlim = c(xmin, xmax), ylim = c(ymin, ymax), ...)
        axis(1)
        axis(2)
        title("Samples not classified to polygons \n")
        # if (moreborders == T) {plot(wrld_simpl, add = T)}
        if (verbose == TRUE) {
            warning("adding polygons")
        }
        if (class(x$polygons) == "list") {
            plota <- function(x) {
                plot(x, add = T, col = "grey60", border = "grey40")
            }
            lapply(x$polygons, plota)
        } else {
            plot(x$polygons, col = "grey60", border = "grey40", add = T, ...)
        }
        if (verbose == TRUE) {
            warning("adding sample points")
        }
        points(dat$XCOOR, dat$YCOOR, cex = 0.5, pch = 3, col = "red", ...)
        box("plot")
    }
}
.NexusOut <- function(dat, verbose = FALSE) {
  if (class(dat) == "list") {
    tablist <- lapply(dat, function(x) x$spec_table)
    for (i in 1:length(tablist)) {
      names(tablist[[i]])[-1] <- paste(names(tablist)[i], names(tablist[[i]][-1]), sep = "_")
    }
    speciestab <- Reduce(function(x, y) merge(x, y, all = TRUE), tablist)
  } else {
    speciestab <- dat$spec_table
  }
  if (verbose == FALSE) {
    sink("species_classification.nex")
  }
  if (verbose == TRUE) {
    sink("species_classification_verbose.nex")
  }
  cat("#NEXUS \n")
  cat("\n")
  cat("begin data; \n")
  cat(paste("\tdimensions ntax=", dim(speciestab)[1], " nchar=", dim(speciestab)[2] - 1, ";", sep = ""))
  cat("\n")
  cat("\tformat datatype=standard symbols=\"01\" gap=-;")
  cat("\n")
  cat("\tCHARSTATELABELS")
  cat("\n")
  if (length(speciestab) == 0) {
    cat("No point fell in any of the polygons specified")
    sink(NULL)
  } else {
    aa <- gsub(" ", "_", names(speciestab)[-1])
    aa <- gsub("&", "_", aa)
    aa <- gsub("__", "_", aa)
    aa <- gsub("__", "_", aa)
    bb <- seq(1, length(aa))
    
    cat(paste("\t", bb[-length(bb)], " ", aa[-length(aa)], ",\n", sep = ""))
    cat("\t", paste(bb[length(bb)], " ", aa[length(aa)], ";\n", sep = ""))
    cat("\n")
    cat("\tmatrix\n")
    
    dd <- as.matrix(speciestab[, -1])
    dd[dd > 0] <- 1
    
    if (dim(dd)[2] > 1) {
      dd <- data.frame(dd)
      dd$x <- apply(dd[, names(dd)], 1, paste, collapse = "")
    } else {
      dd <- data.frame(dd, x = dd)
    }
    ff <- gsub(" ", "_", speciestab[, 1])
    
    if (verbose == F) {
      ee <- paste("\t\t", ff, "\t", dd$x, "\n", sep = "")
      cat(ee)
    }
    if (verbose == T) {
      gg <- vector()
      jj <- speciestab[, -1]
      for (i in 1:dim(jj)[2]) {
        hh <- paste(dd[, i], "[", jj[, i], "]", sep = "")
        gg <- data.frame(cbind(gg, hh))
      }
      gg$x <- apply(gg[, names(gg)], 1, paste, collapse = "")
      ee <- paste("\t\t", ff, "\t", gg$x, "\n", sep = "")
      cat(ee)
    }
    cat("\t;\n")
    cat("end;")
    sink(NULL)
  }
}
.OutBarChartPoly <- function(x, prefix, verbose = FALSE, ...) {
    if (verbose == TRUE) {
        cat("Creating barchart per polygon: barchart_per_polygon.pdf. \n")
    }
    pdf(file = paste(prefix, "barchart_per_polygon.pdf", sep = ""), paper = "special", width = 10.7, height = 7.2, onefile = T)
    .BarChartPoly(x, plotout = T, cex.axis = 0.8, ...)
    dev.off()
}
.OutBarChartSpec <- function(x, prefix, verbose = FALSE, ...) {
    if (verbose == TRUE) {
        cat("Creating barchart per species: barchart_per_species.pdf. \n")
    }
    pdf(file = paste(prefix, "barchart_per_species.pdf", sep = ""), paper = "special", width = 10.7, height = 7.2, onefile = T)
    .BarChartSpec(x, plotout = T, mode = "percent", ...)
    dev.off()
}
.OutHeatCoEx <- function(x, prefix, verbose = FALSE, ...) {
    if (verbose == TRUE) {
        cat("Creating coexistence heatplot: heatplot_coexistence.pdf. \n")
    }
    pdf(file = paste(prefix, "heatplot_coexistence.pdf", sep = ""), paper = "special", width = 10.7, height = 7.2, onefile = T)
    .HeatPlotCoEx(x, ...)
    dev.off()
}
.OutMapAll <- function(x, prefix, areanames = "", verbose = FALSE, ...) {
    if (verbose == TRUE) {
        cat("Creating overview map: map_samples_overview.pdf. \n")
    }
    pdf(file = paste(prefix, "map_samples_overview.pdf", sep = ""), paper = "special", width = 10.7, height = 7.2, onefile = T, 
        ...)
    .MapAll(x, ...)
    .MapUnclassified(x, ...)
    dev.off()
}
.OutMapPerPoly <- function(x, prefix, verbose = FALSE, ...) {
    if (verbose == TRUE) {
        cat("Creating map per polygon: map_samples_per_polygon.pdf. \n")
    }
    pdf(file = paste(prefix, "map_samples_per_polygon.pdf", sep = ""), paper = "special", width = 10.7, height = 7.2, onefile = T)
    .MapPerPoly(x, plotout = T)
    dev.off()
}
.OutMapPerSpecies <- function(x, prefix, verbose = FALSE, ...) {
    if (verbose == TRUE) {
        cat("Creating map per species: map_samples_per_species.pdf. \n")
    }
    pdf(file = paste(prefix, "map_samples_per_species.pdf", sep = ""), paper = "special", width = 10.7, height = 7.2, onefile = T)
    .MapPerSpecies(x, plotout = T, ...)
    dev.off()
}
.OutPlotSpPoly <- function(x, prefix, verbose = FALSE, ...) {
    if (verbose == TRUE) {
        cat("Creating species per polygon barchart: number_of_species_per_polygon.pdf. \n")
    }
    pdf(file = paste(prefix, "number_of_species_per_polygon.pdf", sep = ""), paper = "special", width = 10.7, height = 7.2, 
        onefile = T)
    .PlotSpPoly(x, ...)
    dev.off()
}
.perc <- function(x, y){
  x / rowSums(y) * 100
}
.PipSamp <- function(x, columnname, verbose = FALSE) {
    if (class(x) != "spgeoIN") {
        stop(paste("function is only defined for class \"spgeoIN\".\n", "Use ReadPoints() to produce correct input format", sep = ""))
    }
    occ <- SpatialPoints(x$species_coordinates[, c(1, 2)])
    
    if (class(x$polygons) == "SpatialPolygonsDataFrame") {
        liste <- unique(x$polygons@data[, columnname])
        if (all(!is.na(liste)) == F) {
            # liste <- liste[!is.na(liste)]
            x$polygons@data[, columnname] <- as.character(as.vector(unlist(x$polygons@data[, columnname])))
            x$polygons@data[is.na(x$polygons@data[, columnname]), columnname] <- "unnamed"
            x$polygons@data[, columnname] <- as.factor(as.vector(unlist(x$polygons@data[, columnname])))
            warning("area names contain missing data (#N/A). Renamed to  unnamed. This migh cause problems")
        }
        liste <- unique(x$polygons@data[, columnname])
        bid <- data.frame(x$identifier, as.character(rep(NA, length(x$identifier))), stringsAsFactors = F)
        
        for (i in 1:length(liste)) {
            b <- subset(x$polygons, x$polygons@data[, columnname] == liste[i])
            aaa <- SpatialPolygons(slot(b, "polygons"))
            proj4string(aaa) <- CRS("+proj=longlat +datum=WGS84")
            proj4string(occ) <- CRS("+proj=longlat +datum=WGS84")
            rr <- over(occ, aaa)
            rr[rr > 0] <- as.character(liste[i])
            # if(length(which(rr != 'NA')) != 0){ bid[which(rr != 'NA'),2] <- rr[which(rr != 'NA')] }
            if (length(is.na(rr)) != 0) {
                bid[which(!is.na(rr)), 2] <- rr[which(!is.na(rr))]
            }
        }
        names(bid) <- c("identifier", "homepolygon")
        bid$homepolygon <- as.factor(bid$homepolygon)
        class(bid) <- c("spgeodataframe", "data.frame")
        return(bid)
        
    } else {
        # liste <- levels(x$identifier)
        pp <- x$polygons  #[i]
         proj4string(occ) <- proj4string(pp) <- "+proj=longlat +datum=WGS84"
        if (verbose == TRUE) {
            cat("Performing point in polygon test \n")
        }
        pip <- over(occ, pp)
        if (verbose == TRUE) {
            cat("Done \n")
        }
        pip <- data.frame(x$identifier, pip)
        colnames(pip) <- c("identifier", "homepolygon")
        for (i in 1:length(names(x$polygons))) {
            pip$homepolygon[pip$homepolygon == i] <- names(x$polygons)[i]
        }
        pip$homepolygon <- as.factor(pip$homepolygon)
        class(pip) <- c("spgeodataframe", "data.frame")
        return(pip)
    }
}
.PlotSpPoly <- function(x, ...) {
    if (class(x) == "spgeoOUT") {
        num <- length(names(x$polygon_table))
        dat <- sort(x$polygon_table)
        counter <- num/10
        if (length(x$polygon_table) != 0) {
          if (length(x$polygon_table) == 1){
            barplot(as.matrix(dat[1, ]),
                    ylim = c(0, round((max(dat) + max(c(max(dat)/4, 1))), 0)),
                    ylab = "Number of Species per Polygon", 
                    names.arg = names(x$polygon_table),
                    las = 2, ...)
          }else{
            barplot(as.matrix(dat[1, ]), 
                    ylim = c(0, round((max(dat) + max(c(max(dat)/4, 1))), 0)), 
                    ylab = "Number of Species per Polygon", 
                    las = 2, ...)
          }
          box("plot")
        } else {
            cat("No point in any polygon")
        }
    } else {
        stop("this function is only defined for class \"spgeoOUT\"")
    }
}
.rasterSum <- function(x, ras, type = c("div", "abu")) {
    po <- SpatialPoints(x[, 2:3], CRS("+proj=longlat +datum=WGS84"))
    ras_sub <- rasterize(po, ras, fun = "count")
    if (type == "div") {
        ras_sub[ras_sub >= 1] <- 1
    }
    ras_sub[is.na(ras_sub)] <- 0
    return(ras_sub)
}
.Random.seed <- c(403L, 10L, -1648930018L, -1637691944L, -1360844997L, 244027609L, 250264680L, 1688437142L, 1315422481L, -1590103357L, 
    90966738L, -2034009244L, 420605607L, -2132584451L, -983526924L, 1596694234L, 1916756005L, 801402175L, 473227974L, 454578000L, 
    -1346056333L, -27012719L, -780139584L, -131488354L, 664229737L, 1521537947L, 1326930602L, -927648564L, 2089864591L, 477440453L, 
    -1703916164L, -1989659182L, 1356812621L, -532622329L, 857744366L, 1635232456L, 63407563L, 752033833L, -1699641736L, -1423129690L, 
    1527383489L, 1699562067L, -1240921854L, -1518354380L, 1741420983L, -2102563411L, -773904700L, 202918634L, -1119113899L, 
    -2008816913L, 925958390L, 1198481440L, -721296925L, -235838847L, 407866608L, 1530119054L, 712539385L, -1758280053L, 449214138L, 
    -1451198020L, 1630055167L, 1787579989L, -1003890324L, -1580316222L, -298032611L, 1291696791L, -773876738L, 1587474360L, 
    470336411L, 1406945785L, -1336633016L, -302095242L, 1497182641L, 1917050851L, 1221426354L, 1196036356L, -978888121L, -211509987L, 
    -1802310764L, 1197186874L, -250084475L, 322877535L, 1096616678L, 1285566832L, -556246381L, 1124303665L, 1371496928L, -901917058L, 
    2075974921L, -2065486021L, 22439690L, -389605652L, -614038673L, -1028148059L, 393602652L, 879095346L, 949685549L, 407944167L, 
    -1138410738L, 330956008L, -1431103765L, -2044135543L, -1956199656L, 1886508614L, -1575696735L, 1731365171L, -1016403806L, 
    543598484L, -218135017L, -1120259571L, -305281244L, 896697098L, 1090178677L, 6434255L, 459593430L, -59450752L, 764413123L, 
    1995399777L, -771425584L, -741352402L, 2069131737L, 229666283L, -1727240102L, -214259172L, 2005486367L, -1837580427L, -1805621748L, 
    304889826L, -1454177091L, 261980983L, -1722995746L, -243275496L, 1298406523L, 202572313L, 98460328L, 1753040982L, 86983889L, 
    -1124712317L, -54766830L, 811812516L, 1098510951L, 634022845L, 526435636L, 895958682L, 959778021L, 1080766079L, -967297018L, 
    -181535856L, 894207539L, 147887569L, -1071650944L, -1436419490L, 191537577L, -790425509L, -65058710L, -182392564L, -820370609L, 
    -853147643L, -1739170116L, -620234862L, 446211213L, -1235746105L, 551089070L, 464072200L, -401669621L, -629812503L, -123347528L, 
    1337072230L, -1169725695L, -1338355309L, 667937346L, 1089830388L, -252786697L, -926368787L, -1041034108L, 336595626L, -976035819L, 
    -370341713L, -775485770L, -1461332768L, 1278408739L, 1077831489L, 1479936560L, -1323100082L, -765572935L, 204176587L, -646954630L, 
    1273658492L, 582416575L, 1471499669L, -763032148L, 247220098L, 1634585949L, 1550074967L, 1280172862L, 662768504L, 144189531L, 
    602435513L, -1105285624L, -1017578314L, -2066846735L, 1179883043L, 1253590642L, 1963826372L, 5069447L, 26333789L, 374133844L, 
    760135290L, 97107397L, 386580255L, 1877240486L, 1097387056L, 756182739L, 508824817L, -1573304288L, -1399050562L, 1741179593L, 
    -1092763525L, -1537311158L, -1021220180L, 1672371631L, 593702245L, 1257683740L, -1329422478L, 1214706669L, 1890268199L, 
    -409785266L, -1574940504L, 1977728683L, 1838175008L, -965782068L, 2084854752L, 1538034370L, -782480856L, 842507948L, -692835652L, 
    -628655630L, -2118104784L, -1771234140L, -1632335176L, -1188891814L, 1719432592L, 1081819516L, -1510779628L, 2078191938L, 
    -1777485056L, 1275018172L, 167306064L, -663430878L, 1008647496L, -1324169204L, -1473855396L, 461802706L, -860493824L, -166962348L, 
    -1319540248L, 1561747450L, 1553800912L, -763382484L, -1383411468L, -339109998L, -1832502176L, 70036748L, -1950277920L, 542784866L, 
    506294952L, -1877865908L, -1452267652L, 1272838962L, 867100272L, -628675260L, -1256743976L, -1736621958L, 1879586352L, 175682236L, 
    618276692L, -1244917182L, -502254688L, 454734012L, 1518512592L, -749226654L, -1858905112L, 1701871564L, -92285188L, -144515086L, 
    -920236928L, -5629708L, -1341614968L, 1290601786L, 1695310672L, -246641556L, 1430581396L, 1797118290L, 1062389472L, 909877452L, 
    -975310688L, 1730539906L, 2072891624L, -659447572L, -1035747268L, 1226186866L, -257111952L, 488725028L, 963912248L, -239807974L, 
    1849948368L, 560592444L, 1722406036L, 263365762L, -825036864L, 1930900796L, 1388488336L, 426294306L, -1851807736L, -420105076L, 
    -824448996L, -1073691502L, -1317364608L, -1107411948L, 310391592L, -1484714950L, -1028801328L, 1038126188L, -373344588L, 
    39067666L, 1220337632L, 1633784780L, 1470264800L, -898884574L, 1824331048L, 33734924L, 784211516L, 497669362L, 1972353264L, 
    907487044L, 382304344L, -868992454L, 1914382576L, 478393020L, -560762348L, 595152066L, -477428128L, -85731972L, -161796912L, 
    148855842L, 889325736L, -21316212L, -590698052L, 667793458L, 8231744L, -1100257612L, -267611832L, -2108226886L, 767907600L, 
    778529004L, 1117116500L, 412905234L, 837378464L, -1147170868L, -1781287072L, 2049291458L, -81633752L, 1579780524L, -845876420L, 
    -166030478L, 2146244400L, 269797412L, -133636552L, 790573786L, -576608368L, 665742588L, -151939692L, 153260866L, 1958029824L, 
    -2043934660L, -1981814192L, 1813069986L, 1167496008L, -2070461940L, -1990339876L, -2109028270L, 1326746752L, -1329549484L, 
    1217421288L, -1833018630L, 633356112L, -1375251412L, -570306188L, -272869742L, -1851483552L, 905975948L, 1875831136L, -1494300446L, 
    591135912L, 353558220L, 573295100L, -521391182L, -1106371600L, 438045124L, 904793432L, -742721158L, 311187888L, 1677115836L, 
    1542254164L, -1529988670L, 1221700384L, -1811514436L, -726839216L, 1255302114L, -18856088L, -2065989812L, 930180860L, 1222151922L, 
    1285958912L, 934995060L, -2014024184L, 638180794L, 1007510992L, 460132332L, -1713902700L, 2144571602L, -775832480L, 1772082508L, 
    -1348203744L, 752253570L, 370401768L, -1940134932L, 1354200636L, 1891081842L, -1025035408L, 1799269668L, -1558005192L, 1436619162L, 
    -1484749360L, 1799038396L, -602028396L, 632881794L, -826272320L, 322998076L, -324300144L, 221169442L, 447747592L, -1313193204L, 
    -2016863332L, -264903662L, 1668395904L, -416306412L, -604400344L, 1770745914L, -1244556336L, -1507075476L, -419259340L, 
    186865810L, 2028118368L, -1223881396L, -870986272L, 1853619019L, 1556720060L, 627280362L, -1103917537L, -50492247L, -484535234L, 
    -1368068820L, 127784829L, -845376057L, 904515760L, -552251698L, 1775244283L, 1652233309L, -488748150L, -213320792L, -936347375L, 
    -784887805L, -1171188988L, 483667762L, 115264903L, 2014106577L, 1222847830L, -828374220L, -1213137339L, -1714279793L, 75752456L, 
    -147984346L, -186486573L, 1929314869L, 2006630802L, -1941281440L, 1960627401L, -2069237765L, -440250420L, -1619279270L, 
    -1967079601L, 1974922777L, 643353646L, -946990532L, -1387155539L, 1928095831L, 841305504L, -1216291074L, 1203682763L, 555887565L, 
    -749936902L, 1079235896L, 1526087393L, -1226436205L, -384108108L, 2117413570L, 837999319L, -322073823L, -636966170L, -1722070044L, 
    1930324757L, -392105153L, 1455527640L, -124267786L, -1854378749L, 1784448581L, 1815789154L, 494933456L, -1736350983L, 398741547L, 
    -1841750116L, -155856950L, -1959761729L, 382753289L, 527414046L, -1189823796L, -1591397731L, -427745497L, -1445842608L, 
    1106657262L, 1362600091L, 509872893L, -1323863574L, 1842797192L, -534008783L, -1616355677L, -1044736348L, 984627026L, -2037932377L, 
    -844700943L, 1918216758L, -1584421164L, 608389733L, -1796487057L, 794208936L, 1159034886L, 2080716083L, -171083371L, -1557182990L, 
    -56960448L, -771247191L, 1251421723L, -1843615764L, 1803128826L, 848311343L, 1506706489L, -1477413298L, 545687196L, -1626129395L, 
    568655991L, -767293312L, 1044893406L, 390800811L, 253228717L, -310457574L, -820743208L, -1618071487L, -1810309389L, 941457684L, 
    -447556446L, -388214345L, 1901057409L, 604775302L, -901025340L, -927631115L, 1034819935L, 628757432L, -1282525546L, 22741155L, 
    676915813L, -938770558L, 89789680L, -621024103L, 1628603531L, 722701692L, -619594838L, -1678148641L, -1219759767L, -85353090L, 
    -2087092884L, 1349985725L, 346564231L, -1833278864L, -1779496690L, -371419589L, -1002093283L, -1057927222L, -2072225560L, 
    -1403847727L, 481869251L, 75421380L, -453616270L, -892366777L, 610145297L, 508761494L, -750458636L, 856547973L, 596495311L, 
    795083464L, 973276006L, 1408547219L, 986968821L, 117304018L, -911256288L, 561569033L, -1183389125L, -1080344692L, 856202778L, 
    698263183L, 603447769L, -1354492306L, 1969853948L, 590032749L, -1912121065L, 348627019L)
.SpPerPolH <- function(x, verbose = FALSE) {
    if (verbose == TRUE) {
        cat("Calculating species number per polygon. \n")
    }
    numpoly <- length(names(x)[-1])
    if (numpoly == 0) {
        num_sp_poly <- NULL
    } else {
        pp <- x[, -1]
        pp[pp > 0] <- 1
        if (numpoly > 1) {
            num_sp_poly <- data.frame(t(colSums(pp)))
        } else {
            num_sp_poly <- data.frame(sum(pp))
            names(num_sp_poly) <- names(x)[2]
        }
    }
    return(num_sp_poly)
    if (verbose == TRUE) {
        cat("Done")
    }
}
.SpSumH <- function(x, verbose = FALSE, occ.thresh = occ.thresh) {
  liste <- levels(x$homepolygon)
  if (length(liste) == 0) {
    spec_sum <- NULL
  } else {
    spec_sum <- data.frame(identifier = levels(x$identifier))
    for (i in 1:length(liste)) {
      pp <- subset(x, x$homepolygon == liste[i])
      
      kk <- aggregate(pp$homepolygon, by = list(pp$identifier), length)
      names(kk) <- c("identifier", liste[i])
      spec_sum <- merge(spec_sum, kk, all = T)
    }
    spec_sum[is.na(spec_sum)] <- 0
    if (occ.thresh > 0) {
      filtperc <- data.frame(identifier = spec_sum$identifier, apply(spec_sum[, -1], 2, function(x) .perc(x, y = spec_sum[, -1])))
      
      filtperc <- replace(filtperc, is.na(filtperc), 0)
      for (i in 2:length(names(filtperc))) {
        spec_sum[which(filtperc[, i] < occ.thresh), i] <- 0
      }
    }
  }
  return(spec_sum)
} 
.testcordcap <- function(x, reftab, capthresh) {
  capout <- NA
    if (is.na(as.character(unlist(x["country"]))) | is.na(suppressWarnings(as.numeric(as.character(x["XCOOR"])))) |
          is.na(suppressWarnings(as.numeric(as.character(x["YCOOR"])))) | as.character(unlist(x["country"])) == "") {
        capout <- NA
    } else {
        if (nchar(as.character(unlist(x["country"]))) <= 2) {
            loncap <- suppressWarnings(as.numeric(as.character(x["XCOOR"]))) > (subset(reftab, as.character(reftab$ISO2) == 
                as.character(unlist(x["country"])))$capital_lon - capthresh) & suppressWarnings(as.numeric(as.character(x["XCOOR"]))) < 
                (subset(reftab, as.character(reftab$ISO2) == as.character(unlist(x["country"])))$capital_lon + capthresh)
            latcap <- suppressWarnings(as.numeric(as.character(x["YCOOR"]))) > (subset(reftab, as.character(reftab$ISO2) == 
                as.character(unlist(x["country"])))$capital_lat - capthresh) & suppressWarnings(as.numeric(as.character(x["YCOOR"]))) < 
                (subset(reftab, as.character(reftab$ISO2) == as.character(unlist(x["country"])))$capital_lat + capthresh)
        }
        if (nchar(as.character(unlist(x["country"]))) <= 3 & !nchar(as.character(unlist(x["country"]))) <= 2) {
            loncap <- suppressWarnings(as.numeric(as.character(x["XCOOR"]))) > (subset(reftab, as.character(reftab$ISO3) == 
                as.character(unlist(x["country"])))$capital_lon - capthresh) & suppressWarnings(as.numeric(as.character(x["XCOOR"]))) < 
                (subset(reftab, as.character(reftab$ISO3) == as.character(unlist(x["country"])))$capital_lon + capthresh)
            latcap <- suppressWarnings(as.numeric(as.character(x["YCOOR"]))) > (subset(reftab, as.character(reftab$ISO3) == 
                as.character(unlist(x["country"])))$capital_lat - capthresh) & suppressWarnings(as.numeric(as.character(x["YCOOR"]))) < 
                (subset(reftab, as.character(reftab$ISO3) == as.character(unlist(x["country"])))$capital_lat + capthresh)
        }
        if (nchar(as.character(unlist(x["country"]))) > 3) {
            loncap <- T
            latcap <- T
            warning(paste("found country information for", unlist(x["identifier"]), "with more than 3 letters. Change country information to ISO2 or ISO3", 
                sep = " "))
        }
        ifelse(loncap == T & latcap == T, capout <- FALSE, capout <- TRUE)
    }
    return(capout)
}
.testcordcountr <- function(x, reftab, contthresh) {
  contout <- NA
    if (is.na(as.character(unlist(x["country"]))) | is.na(suppressWarnings(as.numeric(as.character(x["XCOOR"])))) | 
          is.na(suppressWarnings(as.numeric(as.character(x["YCOOR"])))) | as.character(unlist(x["country"])) == "") {
        contout <- NA
    } else {
        if (nchar(as.character(unlist(x["country"]))) <= 2) {
            loncont <- suppressWarnings(as.numeric(as.character(x["XCOOR"]))) > (subset(reftab, as.character(reftab$ISO2) == 
                as.character(unlist(x["country"])))$centroid_lon - contthresh) & suppressWarnings(as.numeric(as.character(x["XCOOR"]))) < 
                (subset(reftab, as.character(reftab$ISO2) == as.character(unlist(x["country"])))$centroid_lon + contthresh)
            latcont <- suppressWarnings(as.numeric(as.character(x["YCOOR"]))) > (subset(reftab, as.character(reftab$ISO2) == 
                as.character(unlist(x["country"])))$centroid_lat - contthresh) & suppressWarnings(as.numeric(as.character(x["YCOOR"]))) < 
                (subset(reftab, as.character(reftab$ISO2) == as.character(unlist(x["country"])))$centroid_lat + contthresh)
        }
        if (nchar(as.character(unlist(x["country"]))) == 3) {
            loncont <- suppressWarnings(as.numeric(as.character(x["XCOOR"]))) > (subset(reftab, as.character(reftab$ISO3) == 
                as.character(unlist(x["country"])))$centroid_lon - contthresh) & suppressWarnings(as.numeric(as.character(x["XCOOR"]))) < 
                (subset(reftab, as.character(reftab$ISO3) == as.character(unlist(x["country"])))$centroid_lon + contthresh)
            latcont <- suppressWarnings(as.numeric(as.character(x["YCOOR"]))) > (subset(reftab, as.character(reftab$ISO3) == 
                as.character(unlist(x["country"])))$centroid_lat - contthresh) & suppressWarnings(as.numeric(as.character(x["YCOOR"]))) < 
                (subset(reftab, as.character(reftab$ISO3) == as.character(unlist(x["country"])))$centroid_lat + contthresh)
        }
        if (nchar(as.character(unlist(x["country"]))) > 3) {
            loncont <- T
            latcont <- T
            warning(paste("found country information for", unlist(x["identifier"]), "with more than 3 letters. Change country information to ISO2 or ISO3", 
                sep = " "))
        }
        ifelse(loncont == T & latcont == T, contout <- FALSE, contout <- TRUE)
    }
    return(contout)
}
.WriteTablesSpGeo <- function(x, prefix = "", verbose = FALSE, ...) {
    if (class(x) == "spgeoOUT") {
        if (verbose == TRUE) {
            cat("Writing sample table: sample_classification_to_polygon.txt. \n")
        }
        write.table(x$sample_table, file = paste(prefix, "sample_classification_to_polygon.txt", sep = ""), row.names = FALSE, sep = "\t", ...)
        if (verbose == TRUE) {
            cat("Writing species occurence table: species_occurences_per_polygon.txt. \n")
        }
        write.table(x$spec_table, file = paste(prefix, "species_occurences_per_polygon.txt", sep = ""), row.names = FALSE, sep = "\t", ...)
        if (verbose == TRUE) {
            cat("Writing species number per polygon table: speciesnumber_per_polygon.txt. \n")
        }
        write.table(x$polygon_table, file = paste(prefix, "speciesnumber_per_polygon.txt", sep = ""), row.names = FALSE, sep = "\t", ...)
        if (verbose == TRUE) {
            cat("Writing table of unclassified samples: unclassified samples.txt. \n")
        }
        write.table(x$not_classified_samples, file = paste(prefix, "unclassified samples.txt", sep = ""), row.names = FALSE, sep = "\t", ...)
        if (verbose == TRUE) {
            cat("Writing coexistence tables: species_coexistence_matrix.txt. \n")
        }
        write.table(x$coexistence_classified, file = paste(prefix, "species_coexistence_matrix.txt", sep = ""), row.names = FALSE, sep = "\t", 
            ...)
    } else {
        stop("this function is only defined for class \"spgeoOUT\"")
    }
} 
