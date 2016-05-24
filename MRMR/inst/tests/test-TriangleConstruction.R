context("TriangleConstruction")

test_that("Cumulative measures", {
  myDF = TestDataFrame()
  myTriangle = newTriangle(TriangleData = myDF
                                 , OriginPeriods = AccidentYear
                                 , DevelopmentLags = Month
                                 , Cumulative = TRUE
                                 , StochasticMeasures = c("Paid")
                                 , StaticMeasures = c("EP")
                                 , Verbose = FALSE)
  
  expect_true(is.Triangle(myTriangle))
  
  df = myTriangle@TriangleData
  expect_true("PriorPaid" %in% colnames(df))
  expect_true("CumulativePaid" %in% colnames(df))
  expect_true("IncrementalPaid" %in% colnames(df))
  expect_false("PriorEP" %in% colnames(df))
})

test_that("Incremental measures", {
  myTriangle = newTriangle(TriangleData = TestDataFrame()
                                 , OriginPeriods = AccidentYear
                                 , DevelopmentLags = Month
                                 , Cumulative = FALSE
                                 , StochasticMeasures = c("Paid")
                                 , StaticMeasures = c("EP")
                                 , Verbose = FALSE)
  
  expect_true(is.Triangle(myTriangle))
  
  df = myTriangle@TriangleData
  expect_true("PriorPaid" %in% colnames(df))
  expect_true("CumulativePaid" %in% colnames(df))
  expect_true("IncrementalPaid" %in% colnames(df))
  expect_false("PriorEP" %in% colnames(df))
})

test_that("Multiple stochastic measures", {
  myTriangle = newTriangle(TriangleData = TestDataFrame()
                           , OriginPeriods = AccidentYear
                           , DevelopmentLags = Month
                           , Cumulative = FALSE
                           , StochasticMeasures = c("Paid", "Reported")
                           , StaticMeasures = c("EP")
                           , Verbose = FALSE)
  
  expect_true(is.Triangle(myTriangle))
  
  df = myTriangle@TriangleData
  expect_true("PriorPaid" %in% colnames(df))
  expect_true("CumulativePaid" %in% colnames(df))
  expect_true("IncrementalPaid" %in% colnames(df))
  expect_false("PriorEP" %in% colnames(df))
  
  expect_true("PriorReported" %in% colnames(df))
  expect_true("CumulativeReported" %in% colnames(df))
  expect_true("IncrementalReported" %in% colnames(df))
})