context("extreme cases")

test_that("Example in which one very small particle and a reasonably sized particle are present",
          {
              expect_true(length(chopthin(c(1e-17,1),2)$weights)==2)
              expect_true(length(chopthin(c(1e-17,0.25,0.25,1),4)$weights)==4)
          }
          )
