context("createMap")

test_that("createMap throws errors", {
  
  expect_error(createMap(NULL, metrics=c("1","3","5")),
               "createMap supports 2 or fewer metrics.")
  
  expect_error(createMap(data.frame(name1=character(0)), metrics=c("name2")),
               "Some of the metrics 'name2' are missing from the data.")
  
  expect_error(createMap(data.frame(name1=character(0)), metrics=c("name2","name1")),
               "Some of the metrics 'name2', 'name1' are missing from the data.")
  
  expect_error(createMap(NULL, location=character(0)),
               "Parameter location is not numeric.")
  
  expect_error(createMap(NULL, location=numeric(0)),
               "Length of parameter location must be 2 or 4.")
  
  expect_error(createMap(NULL, location=c(19)),
               "Length of parameter location must be 2 or 4.")
  
  expect_error(createMap(NULL, location=c(20,21,22)),
               "Length of parameter location must be 2 or 4.")
               
  expect_error(createMap(NULL, location=c(20,21,22,24,25)),
                "Length of parameter location must be 2 or 4.")
})