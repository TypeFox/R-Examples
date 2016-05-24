context("3: Output")

breding_ranges = rgdal::readOGR(system.file(package = "rangeMapper",
     "extdata", "wrens", "vector_combined"), "wrens", verbose = FALSE)[1:10, ]
data(wrens)
d = subset(wrens, select = c('sci_name', 'body_size', 'body_mass', 'clutch_size') )
con = ramp("wrens.sqlite", gridSize = 10, spdf = breding_ranges, biotab = d, ID = "sci_name",
            metadata = rangeTraits(), FUN = "median", overwrite = TRUE)


test_that("rangeMap.fetch does not fetch an empty MAP", {

    rangeMap.save(con, biotab = "biotab", biotrait = "body_mass",
        tableName = "x", FUN = function(x) {NA},
        overwrite = TRUE)

    expect_error ( rangeMap.fetch(con, 'x') )

    rm.rangeMapper(con, tableName = 'MAP_x')

    })


test_that("rangeMap.fetch output", {
	expect_that( rangeMap.fetch(con), is_a('SpatialPixelsRangeMap'))
	expect_that( rangeMap.fetch(con, spatial = FALSE), is_a('rmap.frame'))
	})


test_that("rangeMapExport produces a tiff", {

	m= rangeMap.fetch(con, 'species_richness')
	rastPath = paste( rangeMap.export(con, dir = tempdir()), 'MAP_species_richness.tiff', sep = .Platform$file.sep)

	expect_that(rgdal::readGDAL(rastPath, silent = TRUE), is_a('SpatialGridDataFrame'))

	})
