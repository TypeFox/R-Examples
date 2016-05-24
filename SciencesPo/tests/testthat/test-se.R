context("se expectations")
pdf(NULL) # suppress generating any PDFs
test_that("computes standard error", {
  dat = c(1, 2.3, 2, 3, 4, 8, 12, 43, -1, -4)
  dat %>% se %>% as.numeric %>%
    expect_equal(4.236, tolerance=.005)
  })



