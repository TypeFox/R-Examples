
context("Ant Data")

test_that("We are able to retreive all ant data correctly", {
	# This will also test that georeferencing works correctly
	ad <- aw_data(genus = "acromyrmex")
	ad2 <- aw_data(genus = "acromyrmex", min_elevation =  400)
	expect_true((nrow(ad$data) - nrow(ad2$data)) > 0)

	# # This will also test that georeferencing works correctly
	x1 <- aw_data(min_date = '2014-01-01', max_date = '2014-03-01')
	x2 <- aw_data(min_date = '2014-01-01', max_date = '2014-03-01', georeferenced = TRUE)
	expect_true(x1$count > x2$count)

	expect_is
	expect_is(ad$data, "data.frame")
	expect_error(aw_data())

	acd <- aw_data(genus = "Platythyrea", country = "Madagascar", limit = 20)
	expect_identical(unique(acd$data$country), "Madagascar")
})

context("Data by specimen id works correctly")

test_that("Specimen collections work correctly", {
	data_by_code <- aw_code(catalogNumber="inb0003695883")
	expect_is(data_by_code, "antweb")
	# # BROKEN
	# genus_list <- aw_unique(rank = "genus")
	# expect_is(genus_list, "data.frame")
	# expect_equal(ncol(genus_list), 1)
	# # End broken
	fail <- aw_data(scientific_name = "auberti levithorax")
	expect_is(fail, "NULL")
	fake_code <- aw_code(occurrenceid = "antweb:inb0003695883sdfsdfds") 
	expect_is(fake_code, "NULL")
}) 


test_that("We can correctly retrieve data by coordinates", {
data_by_loc <- aw_coords(coord = "37.76,-122.45", r = 2)
expect_is(data_by_loc$data, "data.frame")
expect_error(aw_coords())
})

context("Combining results")

test_that("we can combine results correctly", {
x1 <- aw_data(genus = "crematogaster", georeferenced = TRUE)
x2 <- aw_data(genus = "crematogaster", georeferenced = TRUE, offset = 1000)
x12 <- aw_cbind(list(x1, x2))
expect_equal(nrow(x1$data), 1000)
expect_equal(nrow(x2$data), 1000)
expect_equal(nrow(x12$data), 2000)
})

context("Photos")

test_that("Photos work correctly", {
	z <- aw_images(since = 5)
	z1 <- aw_images(since = 5, img_type = "d")
	expect_is(z, "data.frame")
	expect_is(z1, "data.frame")
	expect_equal(unique(z1$img_type), "d")
})

test_that("Distinct works correctly", {
	s <- aw_distinct(rank = "genus", country = "Madagascar")
	unique_genera <- length(unique(s$data$genus))
	expect_equal(nrow(s$data), unique_genera)

})

context("Testing the Leaflet maps")

test_that("Leaflet maps and geoJSON work", {
ant_data <- aw_data(genus = "acanthognathus", georeferenced = TRUE)
aw_map(ant_data, dest = ".")
expect_true(file.exists("AntWeb_species_map"))
expect_true(file.exists("temp.geojson"))
unlink("temp.geojson")
unlink("AntWeb_species_map/", recursive = TRUE)
})