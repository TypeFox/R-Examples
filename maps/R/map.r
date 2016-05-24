# R data structure for maps is list(x, y, names)
# names is character vector naming the polygons
# (contiguous groups of non-NA coordinates)

# returns the jth contiguous section of non-NA values in vector x
# should consecutive NA's count as one?
subgroup <- function(x, i) {
  n <- length(x)
  breaks <- which(is.na(x))
  if (length(breaks) == 0) {
    starts <- 1
    ends <- n
  } else {
    starts <- c(1, breaks + 1)
    ends <- c(breaks - 1, n)
  }
## AD: EXTRA: allow negative values in i reverse direction
  pl <- function(j) if(j>0) x[starts[j]:ends[j]] else x[ends[abs(j)]:starts[abs(j)]]
  as.numeric(unlist(lapply(i,function(j) c(NA,pl(j))))[-1])
}

sub.polygon <- function(p, i) {
  lapply(p[c("x", "y")], function(x) subgroup(x, i))
}

# FUTURE:
#sub.gondata <- function(p, i) {
#  x <- p$gondata
#  n <- length(x)
#  breaks <- which(is.na(x))
#  if (length(breaks) == 0) {
#    starts <- 1
#    ends <- n
#  } else {
#    starts <- c(1, breaks + 1)
#    ends <- c(breaks - 1, n)
#  }
#  lengths <- starts - ends + 1
#  i2 <- unlist(lapply(i,function(j) c(NA,x[starts[j]:ends[j]])))[-1]
#  xy <- sub.polygon(p, i2)
#}

# returns a sub-map of the named map corresponding to the given regions
# regions is a vector of regular expressions to match to the names in the map
# regions outside of xlim, ylim may be omitted
map.poly <- function(database, regions = ".", exact = FALSE,
                     xlim = NULL, ylim = NULL, boundary = TRUE,
		     interior = TRUE, fill = FALSE, as.polygon = FALSE, namefield="name") {
  if (!is.character(database)) {
    if (!as.polygon) stop("map objects require as.polygon=TRUE")
    if (inherits(database,"Spatial")){
      if (inherits(database,"SpatialPolygons")) the.map <- SpatialPolygons2map(database, 
                                                               namefield=namefield)
      else if (inherits(database,"SpatialLines")) the.map <- SpatialLines2map(database, 
                                                             namefield=namefield)
      else stop("database not supported.")
    } else the.map <- database
    if (identical(regions,".")) {
      # speed up the common case
      nam = the.map$names
      coord <- the.map[c("x", "y")]
    } else {
      # same as mapname()
      if (exact) {
        i = match(regions, the.map$names)
        if (any(is.na(i))) i = NULL
      } else {
## AD TODO: solve the problem for UK vs Ukrain...
## But here you have a simple set of named polygons, not a map database.
## If it comes from anything but "world", we could break something else.
## We just leave it for now.
      regexp <- paste("(^", regions, ")", sep = "", collapse = "|")
        i <- grep(regexp, the.map$names, ignore.case = TRUE, perl=TRUE)
      }
      if (length(i) == 0) stop("no recognized region names")
      nam <- the.map$names[i]
      coord <- sub.polygon(the.map, i)
    }
    coord$range <- c(range(coord$x, na.rm = TRUE), range(coord$y, na.rm = TRUE))
  } else {
    # turn the polygon numbers into a list of polyline numbers
    gon <- mapname(database, regions, exact)
    n <- length(gon)
    if (n == 0) stop("no recognized region names")
    if (is.null(xlim)) xlim <- c(-1e+30, 1e+30)
    if (is.null(ylim)) ylim <- c(-1e+30, 1e+30)
    # turn the polygon numbers into a list of polyline numbers
    line <- mapgetg(database, gon, as.polygon, xlim, ylim)
    if (length(line$number) == 0)
            stop("nothing to draw: all regions out of bounds")
    # turn the polyline numbers into x and y coordinates
    if (as.polygon) {
      coord <- mapgetl(database, line$number, xlim, ylim, fill) 
      # assemble lines into polygons
      gonsize <- line$size
      keep <- rep(TRUE, length(gonsize))
      coord[c("x", "y")] <- makepoly(coord, gonsize, keep)
    }
    else {
      l <- abs(line$number)
      if (boundary && interior) l <- unique(l)
      else if (boundary) l <- l[!match(l, l[duplicated(l)], FALSE)]
      else l <- l[duplicated(l)]
      coord <- mapgetl(database, l, xlim, ylim, fill)
      if (length(coord) == 0)
              stop("all data out of bounds")
    }
    nam <- line$name
  }
  list(x = coord$x, y = coord$y, range = coord$range, names = nam)
}

.map.range <-
local({
   val <- NULL
   function(new) if(!missing(new)) val <<- new else val
})

map <-
function(database = "world", regions = ".", exact = FALSE,
	 boundary = TRUE, interior = TRUE, projection = "",
	 parameters = NULL, orientation = NULL, fill = FALSE,
	 col = 1, plot = TRUE, add = FALSE, namesonly = FALSE, 
         xlim = NULL, ylim = NULL, wrap = FALSE,
         resolution = if (plot) 1 else 0, type = "l", bg = par("bg"),
         mar = c(4.1, 4.1, par("mar")[3], 0.1), myborder = 0.01, namefield="name", ...)
{
  # parameter checks
  if (resolution>0 && !plot) 
    stop("must have plot=TRUE if resolution is given")
  if (!fill && !boundary && !interior)
    stop("one of boundary and interior must be TRUE")
  doproj <- !missing(projection) || !missing(parameters) || !missing(
          orientation)
  coordtype <- maptype(database)
  if (coordtype == "unknown") 
     stop("missing database or unknown coordinate type")
  if (doproj && coordtype != "spherical") 
    stop(paste(database, "database is not spherical; projections not allowed"))
  # turn the region names into x and y coordinates
  if (is.character(database)) as.polygon = fill
  else as.polygon = TRUE
  coord <- map.poly(database, regions, exact, xlim, ylim, 
                    boundary, interior, fill, as.polygon, namefield=namefield)
  if (is.na(coord$x[1])) stop("first coordinate is NA.  bad map data?")
  if (plot) {
    .map.range(coord$range)
  }
  if (doproj) {
    nam <- coord$names
    coord <- mapproj::mapproject(coord, projection = projection,
			parameters = parameters, orientation = orientation)
    coord$projection = projection
    coord$parameters = parameters
    coord$orientation = orientation
    if (plot && coord$error)
      if (all(is.na(coord$x)))
        stop("projection failed for all data")
      else warning("projection failed for some data")
    coord$names <- nam
  }
  # AD: we do wrapping first: slightly better than when run after the thinning
  #     also now the output data is also wrapped if plot=FALSE
  if (wrap) coord <- map.wrap(coord)
  # do the plotting, if requested
  if (plot) {
    # for new plots, set up the coordinate system;
    # if a projection was done, set the aspect ratio
    # to 1, else set it so that a long-lat square appears
    # square in the middle of the plot
    if (!add) {
      opar = par(bg = bg)
      if (!par("new")) plot.new()
      # xlim, ylim apply before projection
      if (is.null(xlim) || doproj) xrange <- range(coord$x, na.rm = TRUE)
      else xrange <- xlim
      if (is.null(ylim) || doproj) yrange <- range(coord$y, na.rm = TRUE)
      else yrange <- ylim
      if (coordtype != "spherical" || doproj) {
	aspect <- c(1, 1) 
      } else
        aspect <- c(cos((mean(yrange) * pi)/180), 1)
      d <- c(diff(xrange), diff(yrange)) * (1 + 2 * myborder) * aspect
      if (coordtype != "spherical" || doproj) {
        plot.window(xrange, yrange, asp = 1/aspect[1])
      } else {
        # must have par(xpd = FALSE) for limits to have an effect ??!
	par(mar = mar)	# set mai if user-defined mar
	p <- par("fin") -
	  as.vector(matrix(c(0, 1, 1, 0, 0, 1, 1, 0), nrow = 2) %*% par("mai"))
	par(pin = p)
        p <- par("pin")
        p <- d * min(p/d)
        par(pin = p)
        d <- d * myborder + ((p/min(p/d) - d)/2)/aspect
        usr <- c(xrange, yrange) + rep(c(-1, 1), 2) * rep(d, c(2, 2))
        par(usr = usr)
      }
      on.exit(par(opar))
    }
    if (type != "n") {
      # thinning only works correctly if you have polylines from a database
      if (!as.polygon && resolution != 0) {
        pin <- par("pin")
        usr <- par("usr")
        resolution <- resolution * min(diff(usr)[-2]/pin/100)
        coord[c("x", "y")] <- mapthin(coord, resolution)
      }
      if (fill) polygon(coord, col = col, ...)
      else lines(coord, col = col, type = type, ...)
    }
  }
  # return value is names or coords, but not both
  class(coord) = "map"
  value <- if (namesonly) coord$names else coord
  if (plot) invisible(value)
  else value
}

"makepoly" <-
function(xy, gonsize, keep)
{
  # remove NAs and duplicate points so that a set of polylines becomes a set
  # of polygons.
  # xy is a set of polylines, separated by NAs.
  # gonsize is a vector, giving the number of lines to put into each polygon
  # note that a polyline may consist of a single point
  x <- xy$x
  y <- xy$y
  n <- length(x)
  gonsize <- gonsize[ - length(gonsize)]
  discard <- seq(length(x))[is.na(x)]
  if (length(discard) > 0) {
    # locations of (possible) duplicate points
    dups = c(discard - 1, n)
    # only polylines with > 1 point have duplicates
    i = which(diff(c(0, dups)) > 2);
    discard <- c(discard, dups[i]);
  }
  # first part of discard is the NAs, second part is duplicates
  # gonsize tells us which NAs to preserve
  if (length(gonsize) > 0)
    discard <- discard[ - cumsum(gonsize)]
  if (length(discard) > 0) {
    x <- x[ - discard]
    y <- y[ - discard]
  }
  keep <- rep(keep, diff(c(0, seq(length(x))[is.na(x)], length(x))))
  closed.polygon(list(x = x[keep], y = y[keep]))
}
closed.polygon <- function(p) {
  # p is a set of polylines, separated by NAs
  # for each one, the first point is copied to the end, giving a closed polygon
  x = p$x
  y = p$y
  n = length(x)
  breaks <- seq(length(x))[is.na(x)]
  starts <- c(1, breaks + 1)
  ends <- c(breaks - 1, n)
  x[ends + 1] = x[starts]
  y[ends + 1] = y[starts]
  x = insert(x, breaks + 1)
  y = insert(y, breaks + 1)
  list(x = x, y = y)
}
insert <- function(x, i, v = NA) {
  # insert v into an array x, at positions i
  # e.g. insert(1:7, c(2, 5, 8))
  n = length(x)
  new.n = n + length(i)
  m = logical(new.n)
  i.new = i - 1 + seq(length(i))
  m[i.new] = TRUE
  x = x[(1:new.n) - cumsum(m)]
  x[i.new] = v
  x
}
