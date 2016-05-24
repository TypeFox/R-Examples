# ====================
# turn polygons from "maps" into "owin" window for spatstat
# ====================

mapsToOwin <- function(country, database = "worldHires") {

  requireNamespace("mapdata")
  raw <- maps::map(database = database
                   , regions = country
                   , plot = FALSE
                   , fill = TRUE
  )

  cutoff <- which(is.na(raw$x))
  cutoff <- c( 0, cutoff, length(raw$x)+1 )

  coor <- cbind(raw$x, raw$y)

  result <- list()
  for (i in 1:length(raw$names)) {
    result[[i]] <-  coor[(cutoff[i]+1) : (cutoff[i+1]-1), ]
  }
  names(result) <- raw$names

  return( spatstat::owin(poly = result) )
}

# ====================
# turn data from GADM into "owin" windows for spatstat
# ====================

gadmToOwin <- function(country, sub = NULL, level = 0) {

  raw <- raster::getData("GADM", country = country, level = level)

  if (!is.null(sub)) {
    name_objects <- which(grepl("^NAME_", names(raw@data)))
    selection <- sapply(raw@data[name_objects], function(x) {
      which(grepl(sub,x))
    })
    raw <- raw[unlist(selection),]
  }

  return( maptools::as.owin.SpatialPolygons(raw) )

}

# ====================
# turn hull around points into "owin" windows for spatstat
# ====================

hullToOwin <- function(points, shift, alpha) {

  p <- rbind(   points + shift
                , points - shift
                , t(t(points) + c(shift, -shift))
                , t(t(points) + c(-shift, shift))
  )

  hull <- alphahull::ahull(p, alpha =  alpha)

  # turn this into a "owin" window

  hull <- .ah2sp(hull)
  hull <- maptools::as.owin.SpatialPolygons(hull)

  return(hull)
}

# ========================
# code from Andrew Bevan, based on code from Dylan Beaudette
# https://stat.ethz.ch/pipermail/r-sig-geo/2012-March/014409.html
# =========================

.ah2sp <- function(x
                  , increment = 360
                  , rnd = 10
                  , proj4string = sp::CRS(as.character(NA))
                  ) {

  # Extract the edges from the ahull object as a dataframe
  xdf <- as.data.frame(x$arcs)

  # Remove all cases where the coordinates are all the same
  xdf <- subset(xdf, xdf$r > 0)
  res <- NULL

  if (nrow(xdf) > 0){
    # Convert each arc to a line segment
    linesj <- list()
    prevx <- NULL
    prevy <- NULL
    j <- 1

    for (i in 1:nrow(xdf)) {

      rowi <- xdf[i,]
      v <- c(rowi$v.x, rowi$v.y)
      theta <- rowi$theta
      r <- rowi$r
      cc <- c(rowi$c1, rowi$c2)

      # Arcs need to be redefined as strings of points.
      # Work out the number of points to allocate in this arc segment.
      ipoints <- 2 + round(increment * (rowi$theta / 2),0)

      # Calculate coordinates from arc()
      # description for ipoints along the arc.
      angles <- alphahull::anglesArc(v, theta)
      seqang <- seq(angles[1], angles[2], length = ipoints)
      x <- round(cc[1] + r * cos(seqang),rnd)
      y <- round(cc[2] + r * sin(seqang),rnd)

      # Check for line segments that should be joined up
      # and combine their coordinates
      if (is.null(prevx)) {
        prevx <- x
        prevy <- y
      } else if ( x[1] == round(prevx[length(prevx)], rnd)
               && y[1] == round(prevy[length(prevy)], rnd)) {
        if (i == nrow(xdf)) {

          # We have got to the end of the dataset
          prevx <- append(prevx,x[2:ipoints])
          prevy <- append(prevy,y[2:ipoints])
          prevx[length(prevx)] <- prevx[1]
          prevy[length(prevy)] <- prevy[1]
          coordsj <- cbind(prevx,prevy)
          colnames(coordsj) <- NULL

          # Build as Line and then Lines class
          linej <- sp::Line(coordsj)
          linesj[[j]] <- sp::Lines(linej, ID = as.character(j))

        } else {
          prevx <- append(prevx, x[2:ipoints])
          prevy <- append(prevy, y[2:ipoints])
        }

      } else {

        # We have got to the end of a set of lines,
        # and there are several such sets,
        # so convert the whole of this one to a line segment and reset.
        prevx[length(prevx)] <- prevx[1]
        prevy[length(prevy)] <- prevy[1]
        coordsj <- cbind(prevx,prevy)
        colnames(coordsj) <- NULL

        # Build as Line and then Lines class
        linej <- sp::Line(coordsj)
        linesj[[j]] <- sp::Lines(linej, ID = as.character(j))
        j <- j+1
        prevx <- NULL
        prevy <- NULL
      }
    }

    # Promote to SpatialLines
    lspl <- sp::SpatialLines(linesj)

    # Convert lines to polygons
    # Pull out Lines slot and check which lines have start
    # and end points that are the same
    lns <- slot(lspl, "lines")
    polys <- sapply(lns, function(x) {
      crds <- slot(slot(x, "Lines")[[1]], "coords")
      identical(crds[1, ], crds[nrow(crds), ])
    })

    # Select those that do and convert to SpatialPolygons
    polyssl <- lspl[polys]
    list_of_Lines <- slot(polyssl, "lines")
    sppolys <- sp::SpatialPolygons(list(
                sp::Polygons(lapply(list_of_Lines, function(x) {
                  sp::Polygon(slot(slot(x, "Lines")[[1]], "coords"))
                }), ID = "1")), proj4string = proj4string)

    # Create a set of ids in a dataframe
    # then promote to SpatialPolygonsDataFrame
    hid <- sapply(slot(sppolys, "polygons")
                    , function(x) slot(x, "ID"))
    areas <- sapply(slot(sppolys, "polygons")
                    , function(x) slot(x, "area"))
    df <- data.frame(hid,areas)
    names(df) <- c("HID","Area")
    rownames(df) <- df$HID
    res <- sp::SpatialPolygonsDataFrame(sppolys, data=df)
    res <- res[which(res@data$Area > 0),]
  }
  return(res)
}
