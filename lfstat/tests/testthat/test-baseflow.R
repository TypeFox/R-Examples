context("Baseflow Separation")

# Load data set with baseflow computed by Tallaksen
ray <- read.csv2("tallaksen-ray-baseflow.csv")
ray$baseflow <- round(ray$baseflow, 3)
ray$time <- as.Date(ray$time, "%d.%m.%Y")

# compute baseflow with baseflow()
# vector of base flows is as long as input
ray$bf <- baseflow(ray$discharge)


ray96 <- ray[format(ray$time, "%Y") == "1996", c("discharge", "bf")]


ng <- read.csv2("tallaksen-ngaruroro-baseflow.csv")

expect_equal2 <- function(object, expected, tolerance = 1e-10) {
  expect_true(all(abs(object - expected) < tolerance || (is.na(object) & is.na(expected))) )
}

test_that("base flow gets computed correctly", {
  # values according to Tallaksen and van Lanen.
  expect_equal2(ray$baseflow, baseflow(ray$discharge), tolerance = 1e-3)
  expect_equal2(ng$baseflow,  baseflow(ng$discharge),  tolerance = 1e-3)

  # aggregated base flows for river Ray according to Tallaksen and van Lanen.
  # these are mean flow totals per day, not per year as written
  avg <- round(colSums(ray96[, c("discharge", "bf")]), 2)
  expect_equal(as.numeric(avg), c(19.93, 4.03))

  # base flow before first turning point must be NA
  expect_equal(head(ray$baseflow, 15), rep(NA_integer_, 15))

  # base <= discharge
  expect_true(all(ray$baseflow <= ray$discharge, na.rm = T))


})
