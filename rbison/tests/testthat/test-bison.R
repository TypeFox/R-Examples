# tests for bison fxn in taxize
context("bison")

test_that("bison returns the correct value", {
  skip_on_cran()

  out1 <- bison(species="Bison bison", count=5)
  out2 <- bison(species="Canis latrans", type="scientific_name", count=5)
  out3 <- bison(species="Tufted Titmouse", type="common_name")
  out4 <- bison(species = "Phocoenoides dalli dalli", count = 10)
  out5 <- bison(species="Aquila chrysaetos", count=100)
  
  # values
	expect_that(out1$points$name[1], equals("Bison bison"))
	expect_that(out2$points$name[1], equals("Canis latrans"))
	expect_that(out3$points$name[1], equals("Baeolophus bicolor"))
	expect_that(out4$points$name[1], equals("Phocoenoides dalli dalli"))
	expect_that(out5$points$name[1], equals("Aquila chrysaetos"))

	# class
	expect_that(out1, is_a("bison"))
	expect_that(out2, is_a("bison"))
	expect_that(out3, is_a("bison"))
	expect_that(out4, is_a("bison"))
	expect_that(out5, is_a("bison"))
  
	expect_that(out1$summary, is_a("data.frame"))
	expect_that(out1$points, is_a("data.frame"))
	expect_that(out4$counties, equals(NULL))
})
