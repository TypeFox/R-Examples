context("Inputs")



test_that("BEQI-2 data can be read.", {
	expect_that(readBEQI(tempfile()), throws_error())
})



test_that("Redundant spaces will be stripped.", {
	expect_that(stripSpaces("  leading"), is_identical_to("leading"))
	expect_that(stripSpaces("trailing "), is_identical_to("trailing"))
	expect_that(stripSpaces(" Hello   World  "), is_identical_to("Hello World"))
})



test_that("Harmonization removes case-sensitve duplicates.", {
    expect_that(
        harmonize(c("foo", "Foo", "bar", "Foo", "bar")), 
        is_identical_to(c("Foo", "Foo", "bar", "Foo", "bar"))
    )
    expect_that(
        length(unique(harmonize(c("foo", "bar", "Foo", "Bar")))), 
        is_identical_to(2L)
    )
})