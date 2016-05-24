context("test delScattering")
test_that("test EEM class is returned",
          {
              data(applejuice)
              expect_that(delScattering(applejuice), is_a("EEM"))
          })