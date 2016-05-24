context("cut")

test_that("cut works the same by default (with certain sep and format_fun)", {
  x <- runif(10)
  breaks <- seq(0, 1, by = 0.25)
  for (right in c(FALSE, TRUE)) {
    for (include.lowest in c(FALSE, TRUE)) {
      eval(bquote(
        expect_identical(
          cut_format(x, breaks,
                     right = .(right), include.lowest = .(include.lowest),
                     format_fun = formatC, sep = ","),
          cut(x, breaks, right = .(right), include.lowest = .(include.lowest)))
      ))
    }
  }
})

test_that("cut is better", {
  x <- seq(0.1, 0.9, by = 0.2)
  breaks <- seq(0, 1, by = 0.25)
  expect_identical(
    cut_format(x, breaks),
    structure(
      c(1L, 2L, 2L, 3L, 4L),
      .Label = c("(0.00, 0.25]", "(0.25, 0.50]", "(0.50, 0.75]", "(0.75, 1.00]"),
      class = "factor")
  )
})

test_that("custom parentheses", {
  x <- seq(0.1, 0.9, by = 0.2)
  breaks <- seq(0, 1, by = 0.25)
  expect_identical(
    cut_format(x, breaks, paren = c("<", "{", ">", "}")),
    structure(
      c(1L, 2L, 2L, 3L, 4L),
      .Label = c("<0.00, 0.25}", "<0.25, 0.50}", "<0.50, 0.75}", "<0.75, 1.00}"),
      class = "factor")
  )
})

test_that("custom parentheses", {
  expect_error(cut_format(1:10, breaks = 10L), "breaks")
})

test_that("cut_format has same interface", {
  common_names <- intersect(names(formals(cut.default)), names(formals(cut_format)))
  expect_true(all(diff(match(common_names, names(formals(cut.default)))) > 0))
  expect_true(all(diff(match(common_names, names(formals(cut_format)))) > 0))
  expect_identical(setdiff(names(formals(cut.default)), common_names),
                   c("labels", "dig.lab"))
  expect_identical(setdiff(names(formals(cut_format)), common_names),
                   c("format_fun", "sep", "paren"))
})
