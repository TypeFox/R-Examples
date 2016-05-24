
context("2: Process Ranges")

f = system.file(package = "rangeMapper", "extdata", "wrens", "vector_combined")
r = rgdal::readOGR(f, "wrens", verbose = FALSE)[1:10, ]

test_that("reprojecting on the fly", {

	dbcon = rangeMap.start(file = "wrens.sqlite", dir = tempdir(), overwrite = TRUE)
	global.bbox.save(con = dbcon, bbox = f, p4s = CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
	gridSize.save(dbcon, gridSize = 2000000)
	canvas.save(dbcon)

	expect_warning( # because of re-projecting
		processRanges(con = dbcon, spdf = r, ID = "sci_name")  )

	})

test_that("ONE SpPolyDF NO metadata", {

	dbcon = rangeMap.start(file = "wrens.sqlite", dir = tempdir(), overwrite = TRUE)
	global.bbox.save(con = dbcon, bbox = r)
	gridSize.save(dbcon, gridSize = 10)
	canvas.save(dbcon)

	processRanges(con = dbcon, spdf = r, ID = "sci_name")

	expect_true( rangeMap.save(dbcon) )

	expect_that(rangeMap.fetch(dbcon) , is_a("SpatialPixelsRangeMap") )

	expect_that(rangeMap.fetch(dbcon, spatial = FALSE) , is_a("data.table") )


	})

test_that("ONE SpPolyDF WITH metadata", {

	dbcon = rangeMap.start(file = "wrens.sqlite", dir = tempdir(), overwrite = TRUE)

	global.bbox.save(con = dbcon, bbox = r)
	gridSize.save(dbcon, gridSize = 10)
	canvas.save(dbcon)

	metadataFunctions = rangeTraits( simpleRange =  function(x) nrow(vertices(x)) < 10 )

	processRanges(con = dbcon, spdf = r, ID = "sci_name", metadata = metadataFunctions )

	expect_true( rangeMap.save(dbcon) )

	expect_that(rangeMap.fetch(dbcon) , is_a("SpatialPixelsRangeMap") )

	expect_that(rangeMap.fetch(dbcon, spatial = FALSE) , is_a("data.table") )

	expect_more_than(nrow( dbGetQuery(dbcon, 'SELECT * from metadata_ranges') ), 0 )

	})

test_that("MULTIPLE SpPolyDF-s WITH metadata", {

	wd = tempdir()
	rfiles = list.files( system.file(package = "rangeMapper", "extdata", "wrens", "vector") , full.names = TRUE )
	rfiles = rfiles[grep('Odontorchilus', rfiles)]
	file.copy(rfiles, wd )

	dbcon = rangeMap.start(file = "wrens.sqlite", overwrite = TRUE, dir = wd)
	global.bbox.save(con = dbcon, bbox = wd)
	gridSize.save(dbcon, gridSize = 10)
	canvas.save(dbcon)

	processRanges(dir = wd, con = dbcon, metadata = rangeTraits())

	expect_true( rangeMap.save(dbcon) )

	expect_that(rangeMap.fetch(dbcon) , is_a("SpatialPixelsRangeMap") )

	expect_that(rangeMap.fetch(dbcon, spatial = FALSE) , is_a("data.table") )

	expect_more_than(nrow( dbGetQuery(dbcon, 'SELECT * from metadata_ranges') ), 0 )

	})

# !! testthat.R hung when run with R cmd check
# test_that("Parallel", {
#
# 	require(doParallel)
# 	cl = makePSOCKcluster(2)
# 	registerDoParallel(cl)
#
# 	dbcon = rangeMap.start(file = "wrens.sqlite", dir = tempdir(), overwrite = TRUE)
# 	global.bbox.save(con = dbcon, bbox = r)
# 	gridSize.save(dbcon, gridSize = 10)
# 	canvas.save(dbcon)
# 	processRanges(con = dbcon, spdf = r, ID = "sci_name", metadata = rangeTraits() )
#
# 	stopCluster(cl)
# 	registerDoSEQ()
#
# 	})
#
















