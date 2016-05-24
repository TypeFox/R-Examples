context("Truncation of input data")

test_that("climate and tree data are correctly truncated to common timespan", {
  climate <- data.frame(
    year = rep(1950:2009, each = 12),
    month = rep(1:12, 60),
    temp = rep(-10 * cos(seq(0, 2*pi, length.out = 12)), 60),
    prec = rep(seq(100, 220, length.out = 12), 60)
    )
  class(climate) <- c("tcclimate", "data.frame")
  chrono <- data.frame(rnorm(100))
  rownames(chrono) <- 1901:2000
  
  expect_that(truncate_input(chrono, climate,
                             NULL, 1, FALSE)$climate$year,
              equals(rep(1950:2000, each = 12)))
  expect_message(truncate_input(chrono, climate,
                                NULL, -1, FALSE))
})

test_that("climate and tree data are correctly truncated to user supplied specs", {
  climate <- data.frame(
    year = rep(1950:2009, each = 12),
    month = rep(1:12, 60),
    temp = rep(-10 * cos(seq(0, 2*pi, length.out = 12)), 60),
    prec = rep(seq(100, 220, length.out = 12), 60)
    )
  class(climate) <- c("tcclimate", "data.frame")
  chrono <- data.frame(rnorm(100))
  rownames(chrono) <- 1901:2000
  
  expect_that(truncate_input(chrono, climate,
                             c(1955, 1998), 1, FALSE)$climate$year,
              equals(rep(1955:1998, each = 12)))
  
  expect_that(truncate_input(chrono, climate,
                             c(1955, 2005), 1, FALSE),
              throws_error("for start dates in current year"))
  
  expect_that(truncate_input(chrono, climate,
                             c(1955, 2005), -1, FALSE),
              throws_error("for start dates in previous year"))
  
  expect_message(truncate_input(chrono, climate,
                             c(1950, 1998), -1, FALSE))
})
