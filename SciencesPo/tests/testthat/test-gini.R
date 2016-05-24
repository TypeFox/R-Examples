context("Gini-Simpson Index expectations.")
pdf(NULL) # suppress generating any PDFs
test_that("The Gini-Simpson Index", {
  dat <- c(778, 815, 857, 888, 925, 930, 965, 990, 1012);
  dat %>% gini.simpson %>%
    expect_equal(0.888, tolerance=.005)
})



context("Weighted Gini Index expectations.")
pdf(NULL) # suppress generating any PDFs
test_that("The Weighted Gini Index", {
  dat <- c(778, 815, 857, 888, 925, 930, 965, 990, 1012);
  dat %>% gini %>%
    expect_equal(0.0468, tolerance=.005)
})

