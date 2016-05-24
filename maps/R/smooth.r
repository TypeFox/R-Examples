# Functions for smoothing or dis-aggregating data over map regions.

# m is a map object
# z is a named vector
# res is resolution of sampling grid
# span is kernel parameter (larger = smoother)
#   span = Inf is a special case which invokes cubic spline kernel.
#   span is scaled by the map size, and is independent of res.
# result is a frame
smooth.map <- function(m, z, res = 50, span = 1/10, averages = FALSE,
                       type = c("smooth", "interp"), merge = FALSE) {
  #if(is.data.frame(z)) z = as.named.vector(z)
  if(averages) {
    # turn averages into sums
    z = z * area.map(m, names(z), sqmi=FALSE)
  }
  # sampling grid
  xlim <- range(m$x, na.rm = TRUE)
  ylim <- range(m$y, na.rm = TRUE)
  midpoints <- function(start, end, n) {
    inc <- (end - start)/n
    seq(start + inc/2, end - inc/2, len = n)
  }
  # 2*res is an assumption about aspect ratio (usually true)
  if(length(res) == 1) res = c(2*res, res)
  xs <- midpoints(xlim[1], xlim[2], res[1])
  ys <- midpoints(ylim[1], ylim[2], res[2])
  x <- expand.grid(x = xs, y = ys)
  if(FALSE) {
    # add centroids to the sample points
    xc = apply.polygon(m[c("x", "y")], centroid.polygon)
    # convert m into a matrix
    xc <- t(array(unlist(xc), c(2, length(xc))))
    xc = data.frame(x = xc[, 1], y = xc[, 2])
    x = rbind(x, xc)
  }
  radius = sqrt(diff(xlim)^2 + diff(ylim)^2)/2
  lambda = 1/(span*radius)^2
  #cat("lambda = ", lambda, "\n")
  cell.area = diff(xs[1:2])*diff(ys[1:2])

  r <- factor(map.where(m, x))
  if(merge) {
    # merge regions
    # merge[r] is the parent of region r
    # regions with merge[r] = NA are considered absent from the map
    # (no sample points will be placed there)
    # this can be slow on complex maps
    merge <- names(z)
    merge <- merge[match.map(m, merge)]
    names(merge) <- m$names
    levels(r) <- merge[levels(r)]
  }
  # remove points not on the map
  i <- !is.na(r)
  x <- x[i, ]
  r <- r[i]
  xo = x
  if(TRUE) {
    # kludge - drop regions with no samples
    n = table(r)
    bad = (n == 0)
    newlevels = levels(r)
    newlevels[bad] = NA
    levels(r) = newlevels
  }
  # put z in canonical order, and drop values which are not in the map
  z = z[levels(r)]
  # remove regions not named in z, or where z is NA
  bad = is.na(z)
  z = z[!bad]
  newlevels = levels(r)
  newlevels[bad] = NA
  levels(r) = newlevels
  i <- !is.na(r)
  x <- x[i, ]
  r <- r[i]
  # do all regions have sample points?
  n = table(r)
  if(any(n == 0)) stop(paste(paste(names(n)[n == 0], collapse = ", "), "have no sample points"))
  type <- match.arg(type)
  if(FALSE) {
    # code for these is in 315/util.r
    # most time is spent here
    # w <- switch(type,
    #             mass = gp.smoother(x, x, r, lambda),
    #             smooth = kr.smoother(x, x, r, lambda))
    # z = drop(z %*% w)
    # cbind(x, z = z)
  } else {
    if(type == "smooth") {
      z = kernel.smooth(x, z, xo, lambda, r)
    } else {
      z = gp.smooth(x, z, xo, lambda, r)
    }
    z = z/cell.area
    cbind(xo, z = z)
  }
}

gp.smooth <- function(x, z, xo, lambda, r) {
  # predict a function measured at locations x to new locations xo
  krr = kernel.region.region(x, r, lambda)
  white.z = solve(krr, z)
  kernel.smooth(x, white.z, xo, lambda, r, normalize = FALSE)
}

kernel.smooth <- function(x, z, xo, lambda, region = NULL, normalize = TRUE) {
  # predict a function measured at locations x to new locations xo
  if(!is.matrix(x)) dim(x) <- c(length(x), 1)
  if(!is.matrix(xo)) dim(xo) <- c(length(xo), 1)
  n = nrow(x)
  if(is.null(region)) region = 1:n
  if(length(region) < n) stop("region must have same length as x")
  region = as.integer(region)
  if(any(is.na(region))) stop("region has NAs")
  if(max(region) > length(z)) stop("not enough measurements for the number of regions")
  no = nrow(xo)
  if(normalize) {
    # divide by region sizes
    z = as.double(z/as.numeric(table(region)))
  }
  .C("kernel_smooth", PACKAGE="maps",
     as.integer(n), as.integer(ncol(x)),
     as.double(t(x)), z, as.integer(region),
     as.integer(no), as.double(t(xo)), zo = double(no),
     as.double(lambda), as.integer(normalize))$zo
}

kernel.region.region <- function(x, region, lambda) {
  if(!is.matrix(x)) dim(x) <- c(length(x), 1)
  region = as.integer(region)
  nr = max(region)
  krr = .C("kernel_region_region", PACKAGE="maps",
    as.integer(nrow(x)), as.integer(ncol(x)),
    as.double(t(x)),
    region, as.double(lambda), as.integer(nr), krr = double(nr*nr))$krr
  dim(krr) = c(nr, nr)
  krr
}
kernel.region.x <- function(x, region, z, lambda) {
  if(!is.matrix(x)) dim(x) <- c(length(x), 1)
  if(!is.matrix(z)) dim(z) <- c(length(z), 1)
  region = as.integer(region)
  nr = max(region)
  no = nrow(z)
  krx = .C("kernel_region_x", PACKAGE="maps",
    as.integer(nrow(x)), as.integer(ncol(x)),
    as.double(t(x)), region, as.integer(no), as.double(t(z)),
    as.double(lambda), as.integer(nr), krx = double(nr*no))$krx
  dim(krx) = c(nr, no)
  krx
}
