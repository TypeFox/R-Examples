context("importCdf")

## create NetCDF test file
#
#library("RNetCDF")
#nc <- create.nc("tiny.cdf")
#
#dim.def.nc(nc, "scan_number", 2)
#dim.def.nc(nc, "point_number", 10, unlim=TRUE)
#
#var.def.nc(nc, "scan_index", "NC_INT", "scan_number")
#var.def.nc(nc, "point_count", "NC_INT", "scan_number")
#var.def.nc(nc, "scan_acquisition_time", "NC_DOUBLE", "scan_number")
#
#var.def.nc(nc, "mass_values", "NC_DOUBLE", "point_number")
#var.def.nc(nc, "intensity_values", "NC_INT", "point_number")
#
#var.put.nc(nc, "scan_index", c(0,5))
#var.put.nc(nc, "point_count", c(5, 5))
#var.put.nc(nc, "scan_acquisition_time", c(1, 2))
#
#var.put.nc(nc, "mass_values", 1:5, start=1, count=5)
#var.put.nc(nc, "mass_values", 6:10, start=6, count=5)
#var.put.nc(nc, "intensity_values", 11:15, start=1, count=5)
#var.put.nc(nc, "intensity_values", 16:20, start=6, count=5)
#
#close.nc(nc)
#

test_that("importCdf", {
  path <- normalizePath(system.file(file.path("exampledata", "tiny.cdf"),
                                    package="MALDIquantForeign"))

  if (suppressWarnings(require("RNetCDF", quietly=TRUE))) {
    ## suppress warnings to avoid creation of Rplots.pdf
    expect_error(suppressWarnings(MALDIquantForeign:::.importCdf("tmp.tmp")))

    r <- list(createMassSpectrum(mass=1:5, intensity=11:15,
                                 metaData=list(file=path, number=1,
                                               retentionTime=1, scanIndex=0)),
              createMassSpectrum(mass=6:10, intensity=16:20,
                                 metaData=list(file=path, number=2,
                                               retentionTime=2, scanIndex=5)))

    s <- MALDIquantForeign:::.importCdf(path, verbose=FALSE)

    expect_equal(s, import(path, verbose=FALSE))
    expect_equal(s, importCdf(path, verbose=FALSE))
    expect_equal(s, import(path, type="cdf", verbose=FALSE))

    expect_equal(mass(s[[1]]), 1:5)
    expect_equal(intensity(s[[1]]), 11:15)
    expect_equal(basename(metaData(s[[1]])$file), "tiny.cdf")
    expect_equal(s, r)
  } else {
    expect_error(suppressWarnings(MALDIquantForeign:::.importCdf(path)),
                 "install.packages")
  }
})
