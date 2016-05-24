context("test unfold")
test_that("test matrix class is returned",
          {
              data(applejuice)
              expect_that(unfold(applejuice), is_a("matrix"))
          })