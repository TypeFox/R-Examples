context("Test getMeasurementData")

test_that("Testing errors and warnings are given", {
  skip_on_cran()
  
  expect_error(test <- getMeasurementData(),
               regexp = 'Downloading everything will crash your computer')
  expect_error(test <- getMeasurementData(999),
               regexp = 'All measurement types specified are invalid:')
  expect_warning(test <- getMeasurementData(c(12,999)),
                 regexp = 'Some measurement types unknown:')
  expect_error(test <- getMeasurementData('s'),
               regexp = 'All measurement types specified are invalid:')
  expect_warning(test <- getMeasurementData(c('Dispersal Age','s')),
                 regexp = 'Some measurement types unknown:')
  expect_warning(test <- getMeasurementData(measurementType = 12, MSW93Binomial='Pon', silent=TRUE),
                 regexp = 'No data was returned')
  expect_warning(test <- getMeasurementData(measurementType = 12, MSW93Binomial='Pon', locationData = FALSE, silent=TRUE),
                 regexp = 'No data was returned')
  expect_warning(test <- getMeasurementData(measurementType = 12, MSW93Binomial=c('Tom', 'Dick', 'Harry', 'Petaurus breviceps'), silent=TRUE),
                 regexp = 'There were no results returned for the following species: Tom, Dick, Harry')
  expect_warning(test <- getMeasurementData(measurementType = 12, MSW05Binomial='Pon', silent=TRUE),
                 regexp = 'No data was returned')
  expect_warning(test <- getMeasurementData(measurementType = 12, MSW05Binomial=c('Tom', 'Dick', 'Harry', 'Petaurus breviceps'), silent=TRUE),
                 regexp = 'There were no results returned for the following species: Tom, Dick, Harry')
  expect_error(test <- getMeasurementData(measurementType = 12, MSW05Binomial='Tom',MSW93Binomial='Dick'),
               regexp = 'Cannot filter by MSW05Binomial and MSW93Binomial')  
})

test_that("Testing search by measurement type", {
  skip_on_cran()
  
  expect_is(test <- getMeasurementData(measurementType='Dispersal Age',
                                       silent=TRUE), 'data.frame')
  expect_equal(ncol(test), 32)
  expect_is(test2 <- getMeasurementData(c('Growth Data','Dispersal Age'),
                                        silent=TRUE), 'data.frame')
  expect_true(nrow(test) < nrow(test2))
  expect_is(test3 <- getMeasurementData(measurementType = 'Dispersal Age',
                                       silent = TRUE,
                                       locationData = FALSE), 'data.frame')
  expect_equal(ncol(test3), 17)
  
})

test_that("Testing MSW93binomial searches", {
  skip_on_cran()
  
  expect_is(test <- getMeasurementData(measurementType = 12, MSW93Binomial='Petaurus breviceps', silent=TRUE), 'data.frame')
  expect_equal(ncol(test), 32)
  expect_is(test2 <- getMeasurementData(measurementType = 12, MSW93Binomial=c('Petaurus breviceps','Equus zebra'), silent=TRUE), 'data.frame')
  expect_true(nrow(test) < nrow(test2))  
})

test_that("Testing MSW05binomial searches", {
  skip_on_cran()
  
  expect_is(test <- getMeasurementData(measurementType = 12, MSW05Binomial='Petaurus breviceps', silent=TRUE), 'data.frame')
  expect_equal(ncol(test), 32)
  expect_is(test2 <- getMeasurementData(measurementType = 12, MSW05Binomial=c('Petaurus breviceps','Equus zebra'), silent=TRUE), 'data.frame')
  expect_true(nrow(test) < nrow(test2))  
})

test_that("Testing complex search", {
  skip_on_cran()
  
  expect_is(test <- getMeasurementData(measurementType = c(12,16,1), MSW05Binomial=c('Petaurus breviceps','Equus zebra'), silent=TRUE), 'data.frame')
})

test_that("Testing casting", {
  skip_on_cran()
  
  expect_is(test <- getMeasurementData(measurementType = 12, MSW05Binomial='Petaurus breviceps', cast = FALSE, silent=TRUE), 'data.frame')
  expect_equal(ncol(test), 26)
  expect_is(test2 <- getMeasurementData(measurementType = 12, MSW05Binomial='Petaurus breviceps', silent=TRUE), 'data.frame')
  expect_true(ncol(test)<ncol(test2))
})

test_that("Testing location addition", {
  skip_on_cran()
  
  expect_is(test1 <- getMeasurementData(measurementType = 12, locationData = TRUE, silent=TRUE), 'data.frame')
  expect_equal(ncol(test1), 32)
  expect_warning(test2 <- getMeasurementData(measurementType = 12, locationData = TRUE, locationOnly = TRUE, silent=TRUE),
                 regexp = 'rows have been removed as they are') 
  expect_true(nrow(test1) > nrow(test2))
  expect_is(test3 <- getMeasurementData(measurementType = 12, country='Switzerland', silent=TRUE), 'data.frame')
  expect_true(nrow(test2) > nrow(test3))
  expect_is(test4 <- getMeasurementData(measurementType = 12, StudyUnitId=197, silent=TRUE), 'data.frame')  
  expect_true(nrow(test2) > nrow(test4))
  expect_error(test4 <- getMeasurementData(measurementType = 12, StudyUnitId=197, country = 'India', silent=TRUE),
               regexp = 'Cannot use both StudyUnitId and country at the same time')    
})