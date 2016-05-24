context("Pearson's contingency coefficient expectations.\n Friendly (2000), 'Visualizing Categorical Data', SAS Institute Inc., p. 61.")
pdf(NULL) # suppress generating any PDFs
test_that("The Pearson's contingency coefficient", {
  tab = as.table(rbind(c(5, 29, 14, 16), c(15, 54, 14, 10), c(20,  84, 17, 94), c(68, 119, 26, 7)) );
  tab %>% calc.CC %>%
    expect_equal(0.435, tolerance=.005)
})

