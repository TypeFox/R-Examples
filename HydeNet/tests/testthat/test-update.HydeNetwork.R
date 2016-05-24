context("update.HydeNetwork")

carNet <- HydeNetwork(~gear | mpg * am,
                      data = mtcars) 

test_that("update by adding a node",
{
  expect_that(update(carNet, ~ . + cyl | am),
              not(throws_error()))
})
  
test_that("update and lose a parent",
{
  expect_warning( update(carNet, ~ . - am) )
})
    