make.surface.grid <- function (grid.list)
{
  temp <- as.matrix(expand.grid(grid.list))
  dimnames(temp) <- list(NULL, names(grid.list))
  attr(temp, "grid.list") <- grid.list
  temp
}

as.surface <- function (obj, z, order.variables = "xy")
{
  if (is.list(obj)) {
    grid.list <- obj
  }
  if (is.matrix(obj)) {
    grid.list <- attr(obj, "grid.list")
  }
  hold <- parse.grid.list(grid.list, order.variables = "xy")
  c(hold, list(z = matrix(z, ncol = hold$ny, nrow = hold$nx)))
}

parse.grid.list <- function (grid.list, order.variables = "xy")
{
  M <- length(grid.list)
  gcounts <- unlist(lapply(grid.list, FUN = length))
  xy <- (1:M)[gcounts > 1]
  if (length(xy) > 2) {
    stop("only two components of the grid list\ncan have more than one element")
  }
  if (order.variables == "yx") {
    xy <- xy[2:1]
  }
  nx <- gcounts[xy[1]]
  ny <- gcounts[xy[2]]
  x <- grid.list[[xy[1]]]
  y <- grid.list[[xy[2]]]
  xlab <- names(grid.list)[xy[1]]
  ylab <- names(grid.list)[xy[2]]
  xlab <- ifelse(is.null(xlab), "X", xlab)
  ylab <- ifelse(is.null(ylab), "Y", ylab)
  list(x = x, y = y, nx = nx, ny = ny, xlab = xlab, ylab = ylab,
       xy = xy)
}

calcGridFromxy <- function (xyz, nx = 80, ny = 80, xycols = c(1, 2)) {
  grid.list <- as.list(apply(xyz,2,median))
  xrange <- range(xyz[, xycols[1]])
  yrange <- range(xyz[, xycols[2]])
  grid.list[[xycols[1]]] <- seq(xrange[1], xrange[2], , nx)
  grid.list[[xycols[2]]] <- seq(yrange[1], yrange[2], , ny)
  grid.list
}
