context("Cramer's V coefficient for tables expectations.\n Friendly (2000), 'Visualizing Categorical Data', SAS Institute Inc., p. 61.")
pdf(NULL) # suppress generating any PDFs
test_that("The Cramer's V coefficient for tables", {
  tab = as.table(rbind(c(5, 29, 14, 16),
                       c(15, 54, 14, 10),
                       c(20,  84, 17, 94),
                       c(68, 119, 26, 7)));
  tab %>% calc.CV %>%
    expect_equal(0.279, tolerance=.005)
})

