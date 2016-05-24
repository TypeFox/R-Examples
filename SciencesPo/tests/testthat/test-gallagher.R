context("The Gallagher index of LSq index.")
pdf(NULL) # suppress generating any PDFs
test_that("The Gallagher index", {
# 2012 Queensland state elecion
   pvotes= c(49.65, 26.66, 11.5, 7.53, 3.16, 1.47)
   pseats = c(87.64, 7.87, 2.25, 0.00, 2.25, 0.00)
  gallagher(pvotes, pseats) %>%
    expect_equal(31.16, tolerance=.005)
})


