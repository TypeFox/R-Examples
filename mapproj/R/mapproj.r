".Last.projection"<-
local({
    val <- list(projection = "", parameters = NULL, orientation = NULL)
    function(new) if(!missing(new)) val <<- new else val
     })
"mapproject"<-
function(x, y, projection = "", parameters = NULL, orientation = NULL)
{
  # minka: cleaned up handling of defaults
  r <- NULL
  # LY: change test for list so that x$x format not applied to vector
  #if(!is.null(x$x)) {
  if (is.list(x)) {
    r <- x$range[1:2]
    y <- x$y
    x <- x$x
  }
  if(length(x) != length(y))
    stop("lengths of x and y must match")
  if(is.null(r))
    r <- range(x[!is.na(x)])
  new.projection <- (projection != "")
  if(new.projection) {
    if(is.null(orientation)) orientation = c(90, 0, mean(r))
    else if(length(orientation) != 3)
      stop("orientation argument must have 3 elements")
  }
  else {
    if(nchar(.Last.projection()$projection) == 0) {
      #stop("no previous projection")
      return(list(x = x, y = y))
    }
    p <- .Last.projection()
    projection <- p$projection
    if(is.null(parameters)) parameters <- p$parameters
    else if(length(parameters) != length(p$parameters))
      stop(paste("expecting", length(p$parameters),
                 "parameters for", projection, "projection"))

    if(is.null(orientation)) orientation <- p$orientation
    else if(length(orientation) != 3)
      stop("orientation argument must have 3 elements")
  }
  error <- .C("setproj",
              as.character(projection),
              as.double(parameters),
              as.integer(length(parameters)),
              as.double(orientation),
              error = character(1), PACKAGE = "mapproj")$error
  if(error != "")
    stop(error)
  .Last.projection(list(projection = projection,
                        parameters = parameters,
                        orientation = orientation))
  .C("doproj",
     x = as.double(x),
     y = as.double(y),
     as.integer(length(x)),
     range = double(4),
     error = integer(1),
     NAOK = TRUE, PACKAGE = "mapproj")[c("x", "y", "range", "error")]
}

map.grid <-
function(lim, nx = 9, ny = 9, labels = TRUE, pretty = TRUE, cex = 1,
	col = 4, lty = 2, font = 2, ...) {
  # uses map.wrap from maps package
  pretty.range <-
  function(lim, ...) {
    # like pretty but ensures that the range is identical:
    # range(pretty.range(x)) == range(x)
    x = pretty(lim, ...)
    if(abs(x[1]-lim[1]) > abs(x[2]-lim[1])) x = x[-1]
    n = length(x)
    if(abs(x[n]-lim[2]) > abs(x[n-1]-lim[2])) x = x[-n]
    x[1] = lim[1]; x[length(x)] = lim[2]
    x
  }
  auto.format <-
  function(x) {
    # use the minimal number of digits to make x's unique
    # similar to abbrev
    for(digits in 0:6) {
      s = formatC(x, digits = digits, format = "f")
      if(all(duplicated(s) == duplicated(x))) break
    }
    s
  }
  # by default, use limits of last map
  if(missing(lim)) lim = .map.range()
  if(is.list(lim)) {
    # first argument is a map
    lim <- lim$range
  }
  if(lim[2]-lim[1] > 360) {
    lim[2] <- lim[1] + 360
  }
  if(pretty) {
    x <- pretty.range(lim[1:2], n = nx)
    y <- pretty.range(lim[3:4], n = ny)
  } else {
    x <- seq(lim[1], lim[2], len = nx)
    y <- seq(lim[3], lim[4], len = ny)
  }
  p = mapproject(expand.grid(x = c(seq(lim[1], lim[2], len = 100), NA),
  	y = y))
  p = map.wrap(p)
  lines(p,
        col = col, lty = lty, ...)
  lines(mapproject(expand.grid(y = c(seq(lim[3], lim[4],
 	len = 100), NA), x = x)), col = col, lty = lty, ...)
  if(labels) {
    tx <- x[2]
    xinc <- median(diff(x))
    ty <- y[length(y)-2]
    yinc <- median(diff(y))
    text(mapproject(expand.grid(x = x + xinc*0.05,
 	y = ty + yinc*0.5)),
        labels = auto.format(x), cex = cex, adj = c(0, 0), col = col,
	font=font, ...)
    text(mapproject(expand.grid(x = tx + xinc*0.5,
  	y = y + yinc*0.05)),
         labels = auto.format(y), cex = cex, adj = c(0, 0), col = col,
	 font=font, ...)
  }
}
