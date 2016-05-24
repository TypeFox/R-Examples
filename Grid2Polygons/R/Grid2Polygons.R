Grid2Polygons <- function (grd, zcol=1, level=FALSE, at, cuts=20,
                           pretty=FALSE, xlim=NULL, ylim=NULL, ply=NULL) {

  # Additional functions (subroutines)

  # Find polygon nodes:
  #  Input:  s          - matrix; 2-column table giving start- and end-node
  #                       indexes for each segment in a level
  #  Output: poly.nodes - list; vector components giving node indexes for each
  #                       polygon ring. The status of the polygon as a hole or
  #                       an island is taken from the ring direction, with
  #                       clockwise meaning island, and counter-clockwise
  #                       meaning hole.
  FindPolyNodes <- function (s) {

    # Remove duplicate segments
    id <- paste(apply(s, 1, min), apply(s, 1, max), sep="")
    duplicates <- unique(id[duplicated(id)])
    s <- s[!id %in% duplicates, ]

    # Number of segments in level
    m <- nrow(s)

    # Call C program to define polygon rings
    out <- matrix(.C("define_polygons", as.integer(s[, 1]), as.integer(s[, 2]),
                     as.integer(m), "ans"=integer(m * 2L))$ans, nrow=m, ncol=2L)

    # Place returned array into list object
    poly.nodes <- lapply(unique(out[, 2]), function(i) out[out[, 2] == i, 1])

    # Close polygon by joining the first point to the last point
    poly.nodes <- lapply(poly.nodes, function(i) c(i, i[1]))

    return(poly.nodes)
  }


  # Main program

  # Check arguments
  if (!inherits(grd, "SpatialGridDataFrame"))
    stop("Grid object not of class SpatialGridDataFrame")
  if (is.character(zcol) && !(zcol %in% names(grd)))
    stop("Column name not in attribute table")
  if (is.numeric(zcol) && zcol > ncol(slot(grd, "data")))
    stop("Column number outside bounds of attribute table")
  if (!is.null(ply) && !inherits(ply, c("SpatialPolygons",
                                        "SpatialPolygonsDataFrame",
                                        "gpc.poly")))
    stop("Polygon class is incorrect")

  # Crop grid data using limit arguments
  if (!is.null(xlim) | !is.null(ylim)) {
    if (is.null(xlim))
      xlim <- bbox(grd)[1, ]
    if (is.null(ylim))
      ylim <- bbox(grd)[2, ]
    vertices <- matrix(c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1],
                         ylim[1], ylim[1], ylim[2], ylim[2], ylim[1]),
                         nrow=5, ncol=2)
    ply.box <- SpatialPolygons(list(Polygons(list(Polygon(vertices,
                                                          hole=FALSE)), 1)))
    proj4string(ply.box) <- CRS(proj4string(grd))
    grd[[zcol]][is.na(over(grd, ply.box))] <- NA
  }

  # Determine break points
  if (level) {
    if (missing(at)) {
      zlim <- range(grd[[zcol]], finite=TRUE)
      if (pretty)
        at <- pretty(zlim, cuts)
      else
        at <- seq(zlim[1], zlim[2], length.out=cuts)
    }
    zc <- at[1:(length(at) - 1L)] + diff(at) / 2
    z <- zc[findInterval(grd[[zcol]], at, rightmost.closed=TRUE)]
  } else {
    z <- grd[[zcol]]
  }

  # Define nodes and elements
  grd.par <- gridparameters(grd)
  n <- grd.par$cells.dim[1]
  m <- grd.par$cells.dim[2]
  dx <- grd.par$cellsize[1]
  dy <- grd.par$cellsize[2]
  xmin <- grd.par$cellcentre.offset[1] - dx / 2
  ymin <- grd.par$cellcentre.offset[2] - dy / 2
  xmax <- xmin + n * dx
  ymax <- ymin + m * dy
  x <- seq(xmin, xmax, by=dx)
  y <- seq(ymin, ymax, by=dy)
  nnodes <- (m + 1L) * (n + 1L)
  nelems <- m * n
  nodes <- 1L:nnodes
  elems <- 1L:nelems
  coords <- cbind(x=rep(x, m + 1L), y=rep(rev(y), each=n + 1L))
  n1 <- c(sapply(1L:m, function(i) seq(1L, n) + (i - 1L) * (n + 1L)))
  n2 <- n1 + 1L
  n4 <- c(sapply(1L:m, function(i) seq(1L, n) + i * (n + 1L)))
  n3 <- n4 + 1L
  elem.nodes <- cbind(n1, n2, n3, n4)

  # Define segments in each element
  nsegs <- nelems * 4L
  segs <- matrix(data=NA, nrow=nsegs, ncol=4,
                 dimnames=list(1L:nsegs, c("elem", "a", "b", "z")))
  segs[, 1] <- rep(1:nelems, each=4)
  segs[, 2] <- c(t(elem.nodes))
  segs[, 3] <- c(t(elem.nodes[, c(2, 3, 4, 1)]))
  segs[, 4] <- rep(z, each=4)
  segs <- na.omit(segs)

  # Identify levels (or unique values)
  levs <- sort(unique(na.omit(z)))

  # Find polygon nodes for each level
  fun <- function(i) FindPolyNodes(segs[segs[, "z"] == i, c("a", "b")])
  poly.nodes <- lapply(levs, fun)

  # Build lists of Polygon objects
  fun <- function(i) lapply(i, function(j) Polygon(coords[j, ]))
  poly <- lapply(poly.nodes, fun)

  # Build list of Polygons objects
  ids <- make.names(1L:length(poly), unique=TRUE)
  fun <- function(i) Polygons(poly[[i]], ID=ids[i])
  polys <- lapply(1L:length(poly), fun)

  # Convert to SpatialPolygons object, add datum and projection
  sp.polys <- SpatialPolygons(polys, proj4string=CRS(proj4string(grd)))

  # Convert to SpatialPolygonsDataFrame object, add data frame of levels
  d <- data.frame(z=levs, row.names=row.names(sp.polys))
  sp.polys.df <- SpatialPolygonsDataFrame(sp.polys, data=d, match.ID=TRUE)

  # Crop SpatialPolygonsDataFrame object using polygon argument
  if (!is.null(ply)) {
    if (!inherits(ply, "gpc.poly"))
      ply <- as(ply, "gpc.poly")
    is.included <- c()
    ids <- row.names(sp.polys.df)
    p1 <- as(sp.polys.df, "gpc.poly")
    p2 <- list()
    for (i in seq(along=p1)) {
      p <- intersect(p1[[i]], ply)
      is.included[i] <- length(p@pts) > 0
      if (is.included[i]) {
        s.plys <- as(p, "SpatialPolygons")
        l.plys <- lapply(slot(s.plys, "polygons"),
                         function(j) {slot(j, "ID") <- ids[i]; j})
        p2 <- c(p2, l.plys)
      }
    }
    if (length(p2) > 0) {
      p2 <- SpatialPolygons(p2, proj4string=CRS(proj4string(grd)))
      d <- data.frame(z=sp.polys.df[[zcol]][is.included],
                      row.names=ids[is.included])
      sp.polys.df <- SpatialPolygonsDataFrame(p2, data=d)
    }
  }

  return(sp.polys.df)
}
