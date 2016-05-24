map.where <- function(database = "world", x, y, ...)
{
  if(missing(y)) {
    if(is.matrix(x)) { y <- x[, 2]; x <- x[, 1] }
    else if(is.list(x) && !is.null(x$y)) { y <- x$y; x <- x$x }
    else { y <- x[[2]]; x <- x[[1]] }
  }
  if(is.character(database)) {
    dbname <- paste(database, "MapEnv", sep = "")
    # data(list = dbname)
    mapbase <- paste(Sys.getenv(get(dbname)), database, sep = "")
    gon <- .C("map_where", PACKAGE="maps",
       as.character(mapbase),
       as.double(x),
       as.double(y),
       as.integer(length(x)),
       integer(length(x)))[[5]]
    # this must be database, not mapbase
    nam <- mapname(database, ".")
    gon[gon == 0] = NA
    names(nam)[gon]
  }
  else {
    if (inherits(database,"SpatialPolygons")) {
      database <- SpatialPolygons2map(database, ...)
    }
    if(num.polygons(database) != length(database$names))
      stop("map object must have polygons (fill=TRUE)")
    n = length(database$x)
    result = .C("map_in_polygon", PACKAGE="maps",
       as.double(database$x), as.double(database$y), as.integer(n),
       as.double(x), as.double(y), as.integer(length(x)),
       integer(length(x)), NAOK = TRUE)[[7]]
    result[result == 0] = NA
    database$names[result]
  }
}
as.matrix.polygon <- function(x, ...) {
  p = x
  if(is.null(p)) return(p)
  if(is.list(p) && !is.data.frame(p)) p <- cbind(p$x, p$y)
  p
}
in.one.polygon <- function(p, x) {
  # returns a logical vector, whose length is nrow(x)
  if(is.null(p)) return(NA)
  p <- as.matrix.polygon(p)
  if(is.list(x) && !is.data.frame(x)) x <- cbind(x$x, x$y)
  if(is.vector(x)) dim(x) <- c(1, 2)
  # p and x are matrices
  .C("map_in_one_polygon", PACKAGE="maps",
     as.double(p[, 1]), as.double(p[, 2]), as.integer(nrow(p)),
     as.double(x[, 1]), as.double(x[, 2]), as.integer(nrow(x)),
     logical(nrow(x)), as.integer(TRUE))[[7]]
}
in.polygon <- function(p, x) {
  # returns a logical vector, whose length is nrow(x)
  if(is.null(p)) return(NA)
  p <- as.matrix.polygon(p)
  if(is.list(x) && !is.data.frame(x)) x <- cbind(x$x, x$y)
  if(is.vector(x)) dim(x) <- c(1, 2)
  # p and x are matrices
  .C("map_in_polygon", PACKAGE="maps",
     as.double(p[, 1]), as.double(p[, 2]), as.integer(nrow(p)),
     as.double(x[, 1]), as.double(x[, 2]), as.integer(nrow(x)),
     logical(nrow(x)), NAOK = TRUE)[[7]] > 0
}

# polygon is not assumed closed
area.polygon <- function(p) {
  if(is.null(p)) return(NA)
  p <- as.matrix.polygon(p)
  n <- nrow(p)
  x1 <- p[, 1]
  i2 <- c(n, 1:(n - 1))
  x2 <- p[i2, 1]
  y1 <- p[, 2]
  y2 <- p[i2, 2]
  0.5*abs(sum(x1*y2 - x2*y1))
}

centroid.polygon <- function(p) {
  if(is.null(p)) return(c(NA, NA))
  p <- as.matrix.polygon(p)
  n <- nrow(p)
  x1 <- p[, 1]
  i2 <- c(n, 1:(n - 1))
  x2 <- p[i2, 1]
  y1 <- p[, 2]
  y2 <- p[i2, 2]
  a <- x1*y2 - x2*y1
  s <- sum(a)*3
  if(s == 0) c(mean(x1), mean(y1))
  else c(sum((x1 + x2)*a)/s, sum((y1 + y2)*a)/s)
}

# applies fun to all sub-polygons of p
apply.polygon <- function(p, fun, names. = NULL) {
  if(is.null(p)) return(p)
  if(is.null(names.) && !is.null(p$names)) names. = p$names
  p = as.matrix.polygon(p)
  n <- nrow(p)
  breaks <- (1:n)[is.na(p[, 1])]
  starts <- c(1, breaks + 1)
  ends <- c(breaks - 1, n)
  m <- length(starts)
  result <- list()
  for(i in 1:m) {
    this.p = if(ends[i] >= starts[i]) p[starts[i]:ends[i], ] else NULL
    result[[i]] <- fun(this.p)
  }
  names(result) = names.
  result
}

num.polygons <- function(p) {
  if(is.list(p)) 1 + sum(is.na(p$x))
  else 1 + sum(is.na(p[, 1]))
}

# range.polygon <- function(..., na.rm = FALSE) {
#   p <- as.list.polygon(...)
#   lapply(p[c("x", "y")], range, na.rm = na.rm)
# }

map.text <- function(database, regions = ".", exact = FALSE, labels,
    cex = 0.75, add = FALSE, move = FALSE, ...) {
  if(!add) map(database=database, regions=regions, exact=exact, ...)
  # get polygons
  cc = match.call(expand.dots=TRUE)
  cc[[1]] = as.name("map")
  cc$fill = TRUE
  cc$plot = FALSE
  cc$regions = regions
  cc$exact = exact
  cc$move = cc$add = cc$cex = cc$labels = NULL
  cc$resolution = 0
  m = eval(cc)
  if(missing(labels)) {
    labels = gsub(".*,", "", m$names)
  }
  if(num.polygons(m) != length(labels))
    stop("map object must have polygons (fill=TRUE) and equal number of labels")
  x = apply.polygon(m, centroid.polygon)
  # convert m into a matrix
  x <- t(array(unlist(x), c(2, length(x))))
  if(move) {
# AD: this option should probably be removed, as the code is not available
    # require(mining)
    move.collisions2 <- get("move.collisions2")	# to prevent check NOTE
    w = strwidth(labels, units = "inches", cex = cex)
    h = strheight(labels, units = "inches", cex = cex)
    x = move.collisions2(x[, 1], x[, 2], w, h)
  }
  # want to omit map-specific options here (like "exact")
  text(x, labels, cex = cex, ...)
  invisible(m)
}

identify.map <- function(x, n = 1, index = FALSE, ...) {
  # identify polygons in a map
  # must click near the center of the polygon
  m = x
  if(!is.list(m)) stop("must provide a map object")
  if(num.polygons(m) != length(m$names))
    stop("map object must have polygons (fill=TRUE)")
  x = apply.polygon(m, centroid.polygon)
  x <- t(array(unlist(x), c(2, length(x))))
  i = identify(x[, 1], x[, 2], labels = m$names, n = n, ...)
  if(index) i else m$names[i]
}

area.map <- function(m, regions = ".", sqmi=TRUE, ...) {
  # returns the areas of given regions,
  # combining the areas of all regions which match.
  if(!is.list(m)) stop("must provide a map object")
  if(num.polygons(m) != length(m$names))
    stop("map object must have polygons (fill=TRUE)")
  proj = m$projection
  m = map.poly(m,regions,as.polygon=TRUE,...)
  area = unlist(apply.polygon(m, area.polygon))
  merge <- regions[match.map(m, regions, ...)]
  names(merge) <- m$names
  merge = factor(merge, levels = regions)
  area = drop(indicators.factor(merge) %*% area)
  areaSqMiles <- function(proj) {
    # returns a factor f such that f*area.map() is in square miles.
    if(is.null(proj)) proj = "no projection"
    if(proj %in% c("mollweide","azequalarea","aitoff")) 2*15732635
    else if(proj %in% c("sinusoidal","bonne","cylequalarea","albers")) 15732635
    else if(proj == "sp_albers") 15745196
    else {
      warning(paste("sq.mile correction unavailable for",proj))
      1
    }
  }
  if(sqmi) area*areaSqMiles(proj) else area
}
indicators.factor <- function(y) {
  # convert a factor into a matrix of indicators
  # result is level by case
  # works if y contains NAs
  r <- array(0, c(length(levels(y)), length(y)), list(levels(y), names(y)))
  for(lev in levels(y)) r[lev, y == lev] <- 1
  r
}

