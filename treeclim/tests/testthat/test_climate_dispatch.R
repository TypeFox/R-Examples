context("Climate data")

test_that("data.frame input of climate data is processed correctly", {
  df <- data.frame(
    year = rep(1950:2009, each = 12),
    month = rep(1:12, 60),
    temp = rep(-10 * cos(seq(0, 2*pi, length.out = 12)), 60),
    prec = rep(seq(100, 220, length.out = 12), 60)
    )

  dfw <- rbind(df, data.frame(
    year = 2010,
    month = 1,
    temp = -10,
    prec = 100
    ))

  dfp1 <- data.frame(cbind(1950:2009, matrix(df$prec, ncol = 12, byrow
                                             = TRUE)))
  dfp2 <- data.frame(cbind(1950:2009, matrix(df$temp, ncol = 12, byrow
                                             = TRUE)))
  dfpw <- rbind(dfp2, c(2008, rep(1, 12)))
  
  expect_that(as_tcclimate(df), is_a("data.frame"))
  expect_that(as_tcclimate(df)$temp,
              equals(rep(-10 * cos(seq(0, 2*pi, length.out = 12)),
                         60)))
  expect_that(as_tcclimate(dfw), throws_error())

  expect_that(as_tcclimate(dfp1), is_a("data.frame"))
  expect_that(as_tcclimate(dfp2)[,3],
              equals(rep(-10 * cos(seq(0, 2*pi, length.out = 12)),
                         60)))
  expect_that(as_tcclimate(list(dfp1, dfp2)), is_a("data.frame"))
  expect_that(as_tcclimate(dfpw), throws_error())
  expect_that(as_tcclimate(list(dfp1, dfpw)), throws_error())
})
