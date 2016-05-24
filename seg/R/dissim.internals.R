.d.adjust<- function(data, nb) {
  
  if (!is.matrix(nb))
    stop("'nb' must be a matrix", call. = FALSE)
  else if (nrow(nb) != ncol(nb))
    stop("'nb' must be a square matrix", call. = FALSE)
  else if (nrow(nb) != nrow(data))
    stop("nrow(nb) must match nrow(data)", call. = FALSE)
  
  if (sum(nb) != 1)
    warning("the sum of all elements in 'nb' does not equal 1", call. = FALSE)
  
  rowsum <- apply(data, 1, sum)
  removeID <- which(rowsum == 0)
  removeL <- length(removeID)
  if (removeL > 0) {
    warning("remove ", removeL, " rows with no population", call. = FALSE)
    rowsum <- rowsum[-removeID]
    data <- data[-removeID,]
    nb <- nb[-removeID, -removeID]
  }
  
  # Black proportions in census tracts
  z <- data[,1] / rowsum
  # Additional spatial component value
  spstr <- 0
  nbvec <- as.vector(nb)
  INDEX <- which(nbvec != 0)
  for (i in 1:length(INDEX)) {
    rowID <- INDEX[i] %% nrow(nb)
    colID <- INDEX[i] %/% nrow(nb)
    if (rowID == 0)
      rowID <- nrow(nb)
    else
      colID <- colID + 1
    spstr <- spstr + (abs(z[colID] - z[rowID]) * nbvec[INDEX[i]])
  }
  as.vector(spstr)
}

.use.spdep <- function(x, data, p2n.args, n2m.args, verbose) {
  speffect <- NA
  if (require(spdep, quietly = TRUE)) {
    if (verbose) {
      message("library 'spdep' appears to be available")
      message("attempting to calculate Morrill's D(adj)")
    }

    if (missing(p2n.args)) {
      p2n.args <- list(pl = x)
    } else {
      if (is.null(p2n.args$pl))
        p2n.args$pl <- x
    }
    grd.nb <- do.call("poly2nb", p2n.args)
    
    if (missing(n2m.args)) {
      n2m.args <- list(neighbours = grd.nb)
    } else {
      if (is.null(n2m.args$neighbours))
        n2m.args$neighbours <- grd.nb
      if (is.null(n2m.args$style))
        n2m.args$style <- "B"
    }
    grd.nb <- do.call("nb2mat", n2m.args)
    grd.nb <- grd.nb / sum(grd.nb)
    speffect <- .d.adjust(data, grd.nb)
  } else if (verbose) {
    message("failed to load 'spdep'")
  }
  speffect
}

.use.spgrass6 <- function(x, data, wVECT.args, v2n.args, verbose) {
  speffect <- rep(NA, 2)
  if (require(spgrass6, quietly = TRUE) & 
        require(rgdal, quietly = TRUE) & 
        require(spdep, quietly = TRUE)) {
    if (verbose) {
      message("library 'spgrass6' and 'rgdal' appear to be available")
      message("attempting to calculate Wong's D(w) and D(s)")
    }
    
    if (!("SpatialPolygonsDataFrame" %in% is(x)))
      x <- SpatialPolygonsDataFrame(x, as.data.frame(data))

    if (missing(wVECT.args)) {
      wVECT.args <- list(SDF = x, vname = "tmp")
    } else {
      if (is.null(wVECT.args$SDF))
        wVECT.args$SDF <- x
      if (is.null(wVECT.args$vname))
        wVECT.args$vname <- "tmp"
    }
    do.call("writeVECT6", wVECT.args)
    
    if (missing(v2n.args)) {
      v2n.args <- list(vname = wVECT.args$vname)
    } else {
      if (is.null(v2n.args$vname))
      v2n.args$vname <- wVECT.args$vname
    }
    sl <- do.call("vect2neigh", v2n.args)
    sl.mat <- listw2mat(sn2listw(sl))
    sl.mat <- sl.mat / sum(sl.mat)  
    speffect[1] <- .d.adjust(data, sl.mat)
    
    A <- unlist(lapply(slot(x, "polygons"), function(z) slot(z, "area")))
    P <- attr(sl, "total") - attr(sl, "external")
    PAR <- P/A
    maxPAR <- max(PAR)
    
    PAR.mat <- matrix(NA, nrow = length(PAR), ncol = length(PAR))
    for (i in 1:length(PAR)) {
      for (j in 1:length(PAR))
        PAR.mat[i,j] <- ((PAR[i] + PAR[j])/2) / maxPAR
    }
    
    speffect[2] <- .d.adjust(data, sl.mat * PAR.mat)   
  } else if (verbose) {
    message("failed to load 'spgrass6' or 'rgdal'")
  }
  
  speffect
}
