context("Test build_rates") 

test_that("build_rates will fail if input is not of the correct type", {
  
  ## if the start.date or end.date is not Date type, there should be an error
  
  expect_error(build_rates(start.date = "2014-08-22", 
                           end.date = "2014-08-23"))
})

test_that("start date shouldn't be later than end date", {
  expect_error(build_rates(start.date = as.Date("2014-08-22"), 
                           end.date = as.Date("2014-08-21")))
})