context("test fold")
test_that("test EEM class is returned",
          {
              data(applejuice)
              unfold_data <- unfold(applejuice)
              expect_that(fold(unfold_data), is_a("EEM"))
          })