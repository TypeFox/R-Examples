library(ncdf4)
library(PCICt)

ncdf4.helpers.test.nc.get.time.series <- function() {
  f <- list(t=nc_open("test1.nc"), not=nc_open("test1.nc", readunlim=FALSE))
  correct.data.ts.test1 <- structure(c(599227200, 599313600, 599400000, 599486400, 599572800),
                                     .Dim = 5L, cal = "365", months = c(31, 28, 31, 30, 31, 30,
                                                               31, 31, 30, 31, 30, 31),
                                     class = "PCICt", dpy = 365, tzone = "GMT", units = "secs")
  correct.data.ts.test2 <- structure(c(599227200, 599313600, 599400000, 599486400, 599572800),
                                     .Dim = 5L, cal = "365", months = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31),
                                     class = "PCICt", dpy = 365, tzone = "GMT", units = "secs",
                                     bounds = structure(c(599184000, 599270400, 599270400, 599356800, 599356800, 599443200, 599443200,
                                       599529600, 599529600, 599616000), .Dim = c(2L, 5L), cal = "365",
                                       months = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31), class = "PCICt", dpy = 365, tzone = "GMT", units = "secs"))
  
  lapply(f, function(f1) {
    checkEquals(nc.get.time.series(f1), correct.data.ts.test1)
    checkEquals(nc.get.time.series(f1, return.bounds=TRUE), correct.data.ts.test2)
    checkEquals(nc.get.time.series(f1, "tasmax"), correct.data.ts.test1)
    checkEquals(nc.get.time.series(f1, time.dim.name="time"), correct.data.ts.test1)
    
    checkException(nc.get.time.series(f1, "foo"))
    checkException(nc.get.time.series(f1, time.dim.name="foo"))
    nc_close(f1)
  })
}

ncdf4.helpers.test.nc.get.proj4.string <- function() {
  f <- nc_open("test1.nc", readunlim=FALSE)
  proj4.string <- "+proj=ob_tran +o_proj=longlat +lon_0=-97 +o_lat_p=42.5 +a=1 +to_meter=0.0174532925199 +no_defs"
  checkEquals(nc.get.proj4.string(f, "tasmax"), proj4.string)
  lat.dat <- ncvar_get(f, "lat")
  lon.dat <- ncvar_get(f, "lon")

  indices <- matrix(c(1, 155, 155, 1, 1, 1, 130, 130), nrow=4, ncol=2)
  colnames(indices) <- c("x", "y")
  
  projected.data <- list(x=f$dim$rlon$vals[indices[,"x"]], y=f$dim$rlat$vals[indices[,"y"]])
  
  latlon.data <- project(projected.data, proj4.string, ellps.default=NA, inverse=TRUE)
  checkEquals((latlon.data$x + 360) %% 360, lon.dat[indices], tolerance=1e-5)
  checkEquals(latlon.data$y, lat.dat[indices], tolerance=1e-5)
  nc_close(f)
}

ncdf4.helpers.test.nc.get.compress.dims <- function() {
  f <- nc_open("test1.nc", readunlim=FALSE)
  compress.dims <- nc.get.compress.dims(f, "tasmax")
  checkEquals(compress.dims, list())
  nc_close(f)                                                                                                                                                                                                                                                                 
}

ncdf4.helpers.test.nc.get.dim.axes <- function() {
  f <- nc_open("test1.nc", readunlim=FALSE)
  dim.axes <- structure(c("Y", "X", "T", NA), .Names = c("rlat", "rlon", "time", "bnds"))
  dim.axes.var <- structure(c("X", "Y", "T"), .Names = c("rlon", "rlat", "time"))
  checkEquals(nc.get.dim.axes(f), dim.axes)
  checkEquals(nc.get.dim.axes(f, "tasmax"), dim.axes.var)
  nc_close(f)
}

ncdf4.helpers.test.nc.get.coordinate.axes <- function() {
  f <- nc_open("test1.nc", readunlim=FALSE)
  coord.axes <- structure(c("X", "Y"), .Names = c("lon", "lat"))
  checkEquals(nc.get.coordinate.axes(f, "tasmax"), coord.axes)
  nc_close(f)
}

ncdf4.helpers.test.nc.get.dim.axes.from.names <- function() {
  f <- nc_open("test1.nc", readunlim=FALSE)
  dim.axes <- structure(c(NA, NA, "T", NA), .Names = c("rlat", "rlon", "time", "bnds"))
  dim.axes.var <- structure(c(NA, NA, "T"), .Names = c("rlon", "rlat", "time"))
  checkEquals(nc.get.dim.axes.from.names(f), dim.axes)
  checkEquals(nc.get.dim.axes.from.names(f, "tasmax"), dim.axes.var)
  nc_close(f)
}

ncdf4.helpers.test.nc.get.dim.names <- function() {
  f <- nc_open("test1.nc", readunlim=FALSE)
  dim.names <- c("rlat", "rlon", "time", "bnds")
  dim.names.var <- c("rlon", "rlat", "time")
  checkEquals(nc.get.dim.names(f), dim.names)
  checkEquals(nc.get.dim.names(f, "tasmax"), dim.names.var)
  nc_close(f)
}

ncdf4.helpers.test.nc.get.variable.list <- function() {
  f <- nc_open("test1.nc", readunlim=FALSE)
  checkEquals(nc.get.variable.list(f), "tasmax")
  checkEquals(nc.get.variable.list(f, min.dims=4), character(0))
  nc_close(f)
}

ncdf4.helpers.test.nc.get.climatology.bounds.var.list <- function() {
  f <- nc_open("test1.nc", readunlim=FALSE)
  checkEquals(nc.get.climatology.bounds.var.list(f), logical(0))
  nc_close(f)
}

ncdf4.helpers.test.nc.get.dim.bounds.var.list <- function() {
  f <- nc_open("test1.nc", readunlim=FALSE)
  checkEquals(nc.get.dim.bounds.var.list(f), structure("time_bnds", .Names = "time"))
  nc_close(f)
}

ncdf4.helpers.test.nc.get.dim.for.axis <- function() {
  f <- nc_open("test1.nc", readunlim=FALSE)
  checkEquals(nc.get.dim.for.axis(f, "tasmax", "X")$name, "rlon")
  nc_close(f)
}

ncdf4.helpers.test.nc.get.var.subset.by.axes <- function() {
  f <- nc_open("test1.nc", readunlim=FALSE)
  all.data <- ncvar_get(f, "tasmax")
  dimnames(all.data) <- list(NULL, NULL, NULL)
  checkEquals(nc.get.var.subset.by.axes(f, "tasmax", list(X=1:4, Y=c(1, 3, 5))), all.data[1:4, c(1, 3, 5), ])
  nc_close(f)
}

## Tests nc.put.var.subset.by.axes, nc.conform.data, nc.copy.atts
ncdf4.helpers.test.write.functions <- function() {
  filename <- tempfile()
  f.in <- nc_open("test1.nc")
  dat <- nc.get.var.subset.by.axes(f.in, "tasmax", list(X=1:4, Y=c(1, 3, 5)))
  x.dim <- ncdim_def("rlon", "degrees", vals=f.in$dim$rlon$vals[c(3:4, 1:2)])
  y.dim <- ncdim_def("rlat", "degrees", vals=f.in$dim$rlat$vals[c(1, 3, 5)])
  t.dim <- ncdim_def("time", "days since 1949-12-01", vals=f.in$dim$time$vals)
  var.list <- list(tasmax=ncvar_def("tasmax", "K", list(t.dim, x.dim, y.dim), 1e20, longname="Daily Maximum Near-Surface Air Temperature"))
  f.out <- nc_create(filename, var.list)
  nc.copy.atts(f.in, "rlat", f.out, "rlat")
  nc.copy.atts(f.in, "rlon", f.out, "rlon")
  nc.copy.atts(f.in, "tasmax", f.out, "tasmax")
  dat.permuted <- nc.conform.data(f.in, f.out, "tasmax", "tasmax", dat, allow.dim.subsets=TRUE)
  nc.put.var.subset.by.axes(f.out, "tasmax", dat.permuted, list())
  nc_sync(f.out)
  nc_close(f.out)
  f.out <- nc_open(filename)

  dat.out <- nc.get.var.subset.by.axes(f.out, "tasmax", list())
  checkEquals(as.numeric(dat.out), as.numeric(dat.permuted))
  checkEquals(dat[1, 2, 3], dat.out[3, 3, 2])

  nc_close(f.in)
  nc_close(f.out)

  unlink(filename)
}
