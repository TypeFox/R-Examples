context("Period of Analysis tests")

test_that("pa", {
  
  year <- 2000
  paStart <- 10
  paLong <- 12
  vectorYear <- c(seq(1999,2001,0.05))
  paIndexWaterYear <- paVector(year, paStart, paLong, vectorYear)
  requestedYears <- vectorYear[paIndexWaterYear]
  expect_equal(requestedYears, c(seq(1999.75,2000.7, 0.05)))
  
  paStart <- 11
  paLong <- 3
  paIndexWinter <- paVector(year, paStart, paLong, vectorYear)
  requestedWinterYears <- vectorYear[paIndexWinter]
  expect_equal(requestedWinterYears, c(seq(1999.85,2000.05, 0.05)))
  
  paStart <- 6
  paLong <- 3
  paIndexSummer <- paVector(year, paStart, paLong, vectorYear)
  requestedSummerYears <- vectorYear[paIndexSummer]
  expect_equal(requestedSummerYears, c(seq(2000.45,2000.65, 0.05)))
  
  paStart <- 10
  paLong <- 3
  paIndexLate <- paVector(year, paStart, paLong, vectorYear)
  endOfYear <- vectorYear[paIndexLate]
  expect_equal(endOfYear, c(seq(2000.75,2000.95, 0.05)))
    
})
