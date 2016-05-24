context("test cutEEM")
test_that("test EEM class is returned",
          {
              data(applejuice)
              expect_that(cutEEM(applejuice), is_a("EEM"))
          })
test_that("test cutEX cannot accept value that does not cover neither min(EX) or max(EX)",
          {
              data(applejuice)
              expect_that(cutEEM(applejuice, cutEX = 300:400), 
                          throws_error("Cannot cut through the middle."))
          })
test_that("test cutEM cannot accept value that does not cover neither min(EM) or max(EM)",
          {
              data(applejuice)
              expect_that(cutEEM(applejuice, cutEM = 400:500), 
                          throws_error("Cannot cut through the middle."))
          })