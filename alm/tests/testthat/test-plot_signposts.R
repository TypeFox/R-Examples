# tests for plot_signposts fxn in alm
context("plot_signposts")

# Plot data from a single identifier gives a bar chart
dat <- alm_signposts(doi="10.1371/journal.pone.0029797")
p <- suppressMessages(plot_signposts(input=dat))

# Plot data from many identifiers gives a line chart
dois <- c('10.1371/journal.pone.0001543','10.1371/journal.pone.0040117',
          '10.1371/journal.pone.0029797','10.1371/journal.pone.0039395')
dat <- alm_signposts(doi=dois)
q <- suppressMessages(plot_signposts(input=dat))

test_that("plot_signposts returns the correct class", {
  expect_that(p, is_a("gg"))
  expect_that(p, is_a("ggplot"))
  expect_that(p$data, is_a("data.frame"))
  expect_that(p$layers[[1]], is_a("proto"))
  expect_that(p$theme, is_a("theme"))
  expect_that(p$coordinates, is_a("cartesian"))
  expect_that(q, is_a("gg"))
  expect_that(q, is_a("ggplot"))
})

test_that("plot_signposts returns the correct dimensions", {
  expect_that(nrow(p$data), equals(4))
  expect_that(nrow(q$data), equals(16))
})
