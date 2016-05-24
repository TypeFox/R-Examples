
context("annotate")

test_that("1 `genocov.files` & subset of SNPs in 1 `snplists.files`", {
  snps <- c("rs2731672", "rs10081087", "rs2287694", "rs7712944")
  
  tab1 <- annotateSNPs(snps)
  tab2 <- annotateSNPs(snps, query.size = 3)
  tab3 <- annotateSNPs(snps, query.size = 2)
  
  expect_true(all(snps %in% tab1$Query))
  expect_true(all(snps %in% tab2$Query))
  expect_true(all(snps %in% tab3$Query))  
})
