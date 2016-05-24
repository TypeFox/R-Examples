# distogram.R
#
# Aron Eklund
#
# A part of the "squash" R package


diamond <- function(x, y = NULL, radius, ...) {
  xy <- xy.coords(x, y)
  xL <- xy$x - radius
  xC <- xy$x
  xR <- xy$x + radius
  yB <- xy$y - radius
  yC <- xy$y
  yT <- xy$y + radius
  n <- length(xL)
  x2 <- rbind(rep.int(NA, n), xC, xL, xC, xR)[-1]
  y2 <- rbind(rep.int(NA, n), yB, yC, yT, yC)[-1]
  polygon(x2, y2, ...)
}        


trianglegram <- function(x, labels = rownames(x), 
    lower = TRUE, diag = FALSE, right = FALSE, 
    add = FALSE, xpos = 0, ypos = 0, xlim, ylim, ...) {
  if(nrow(x) != ncol(x))
    stop("x must be a square matrix")
  n <- nrow(x)
  if(lower) {
    wh <- lower.tri(x, diag = diag)
  } else {
    wh <- upper.tri(x, diag = diag)
  }
  ## x1, y1 = unrotated coordinates
  x1 <- col(x)[wh]
  y1 <- row(x)[wh]
  ## rotated coordinates
  if(right) {
    x2 <- (  y1 - x1) / 2
    y2 <- (x1 + y1) / 2 
  } else {     # right
    x2 <- (- y1 + x1) / 2
    y2 <- (x1 + y1) / 2 
  }
  x2 <- x2 + xpos
  y2 <- y2 + ypos
  if(is.null(labels)) labels <- 1:n
  if(!add) {
    if(missing(xlim)) xlim <- c(-n/2, n/2) + xpos
    if(missing(ylim)) ylim <- c(0.5, n + 0.5) + ypos
    plot(xlim, ylim, type = 'n', 
      axes = FALSE, xlab = '', ylab = '', ...)
  }
  diamond(x2, y2, radius = 0.5, col = x[wh])
  if(diag) offset <- 0.5 else offset <- 0
  if(right) {
    text(0, 1:n, labels = labels, pos = 2, offset = offset, xpd = NA)
  } else {
    text(0, 1:n, labels = labels, pos = 4, offset = offset, xpd = NA)
  }
}  


distogram <- function(x, map, n = 10, base = NA, colFn = heat,
    key = TRUE, title = NA, ...) {
  if(class(x) == 'dist') x <- as.matrix(x)
  stopifnot(nrow(x) == ncol(x))
  if (missing(map)) {
    map <- makecmap(x, n = n, base = base, colFn = colFn)
  }
  trianglegram(cmap(x, map = map), ...)
  if(key) hkey(map = map, title = title)
  invisible(map)
}  


corrogram <- function(...) {
  map <- makecmap(c(-1,1), n = 20, colFn = blueorange, include.lowest = TRUE)
  distogram(..., map = map)
}  


