# tests for alm fxn in alm
context("Testing main functions")

test_that("Main functions return data.frames", {
cc <- of_country_codes()
sp <- of_species_codes()
# landings by country
landings <- of_landings(country = 'CAN')
# landings by species
lbysp <- of_landings(species = "SKJ")
    expect_that(cc, is_a("data.frame"))
    expect_that(sp, is_a("data.frame"))
    expect_that(landings, is_a("data.frame"))
    expect_that(lbysp, is_a("data.frame"))
    expect_equal(ncol(sp), 5)
    expect_equal(ncol(cc), 2)
    expect_equal(ncol(landings), 3)

})


test_that("Ensure that functions fail when presented with bad arguments", {
	expect_error(of_landings(species = "foo"), NULL)
    expect_error(of_landings(country = "foo"), NULL)
    expect_error(of_landings(species = "foo", country = "foo"))
})


test_that("Visualizations are of the right class", {
	test_plot <- fish_plot(of_landings(species = "COD"))
	expect_that(test_plot, is_a("ggplot"))
})
