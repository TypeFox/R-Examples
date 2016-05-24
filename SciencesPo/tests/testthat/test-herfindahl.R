context("Herfindahl Index expectations.")
pdf(NULL) # suppress generating any PDFs
test_that("The Herfindahl Index of concentration", {
  dat <- c(778, 815, 857, 888, 925, 930, 965, 990, 1012);
  dat %>% herfindahl(parameter=1) %>%
    expect_equal(0.1119, tolerance=.005)
})
