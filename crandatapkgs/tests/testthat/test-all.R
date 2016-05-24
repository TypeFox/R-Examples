context("Test return object")

options(repos = c("https://cloud.r-project.org"))
a <- available_data()

test_that("data.frame output", {
    expect_true(is.data.frame(a))
})

test_that("has length", {
    expect_true(nrow(a) > 0)
})

test_that("contains correct packages", {
    expect_true("ggplot2movies" %in% a$Package)
})

test_that("contains correct packages", {
    expect_true(!"MTurkR" %in% a$Package)
})
