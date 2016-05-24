context("Cramer's Phi coefficient for tables expectations.\n Friendly (2000), 'Visualizing Categorical Data', SAS Institute Inc., p. 63.")
pdf(NULL) # suppress generating any PDFs
test_that("The Phi coefficient for tables", {
  tab = as.table(rbind(c(1198, 557), c(1493, 1278)));
  tab %>% calc.Phi %>%
    expect_equal(.143, tolerance=.005)
})




