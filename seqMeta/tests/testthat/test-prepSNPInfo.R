context("prepSNPInfo")


test_that("correct SNPInfo returns as expected", {
  si <- data.frame(Name=c("SNP1", "SNP2"), 
                   gene=c("gene1", "gene1"),
                   wt1=c(.5, .5), 
                   wt2=c(.25, .25),
                   stringsAsFactors=FALSE)
  
  expect_equal(prepSNPInfo(si, "Name", "gene"), si[ , c(2:1)])
  expect_equal(prepSNPInfo(si, "Name", "gene", wt1="wt1"), si[ , c(2:1,3)])
  expect_equal(prepSNPInfo(si, "Name", "gene", wt2="wt2"), si[ , c(2:1,4)])
  expect_equal(prepSNPInfo(si, "Name", "gene", wt1="wt1", wt2="wt2"), si[ , c(2:1,3:4)])
})


test_that("Non-unique weights for a snp but in a different gene", {
  si <- data.frame(Name=c("SNP1", "SNP1"), 
                   gene=c("gene1", "gene2"), 
                   wt1=c(.5, .25), wt2=c(.25, .5), 
                   stringsAsFactors=FALSE)
  
  expect_equal(prepSNPInfo(si, "Name", "gene"), si[ , c(2:1)])
  expect_equal(prepSNPInfo(si, "Name", "gene", wt1="wt1"), si[ , c(2:1,3)])
  expect_equal(prepSNPInfo(si, "Name", "gene", wt2="wt2"), si[ , c(2:1,4)])
  expect_equal(prepSNPInfo(si, "Name", "gene", wt1="wt1", wt2="wt2"), si[ , c(2:1,3:4)])
})


test_that("missing  weight generates error", {
  si <- data.frame(Name=c("SNP1", "SNP2"), 
                   gene=c("gene1", "gene1"), 
                   wt1=c(.5, NA), 
                   wt2=c(.25, .25), 
                   stringsAsFactors=FALSE)
  
  expect_equal(prepSNPInfo(si, "Name", "gene"), si[ , c(2:1)])
  expect_equal(prepSNPInfo(si, "Name", "gene", wt2="wt2"), si[ , c(2:1,4)])
  expect_error(prepSNPInfo(si, "Name", "gene", wt1="wt1"), regexp="missing")
  expect_error(prepSNPInfo(si, "Name", "gene", wt2="wt1"), regexp="missing")
  expect_error(prepSNPInfo(si, "Name", "gene", wt1="wt1", wt2="wt2"), regexp="missing")
  expect_error(prepSNPInfo(si, "Name", "gene", wt1="wt2", wt2="wt1"), regexp="missing")
})


test_that("Non-unique weights for the same snp-gene pair", {
  si <- data.frame(Name=c("SNP1", "SNP1"), 
                   gene=c("gene1", "gene1"), 
                   wt1=c(.5, .25), 
                   wt2=c(.25, .5), 
                   stringsAsFactors=FALSE)
  
  expect_error(prepSNPInfo(si, "Name", "gene", wt1="wt1"), regexp="unique")
  expect_error(prepSNPInfo(si, "Name", "gene", wt2="wt1"), regexp="unique")
  expect_error(prepSNPInfo(si, "Name", "gene", wt1="wt1", wt2="wt2"), regexp="unique")
  expect_error(prepSNPInfo(si, "Name", "gene", wt1="wt2", wt2="wt1"), regexp="unique")
})


test_that("duplicated snps are dropped", {
  si <- data.frame(Name=c("SNP1", "SNP1"), 
                   gene=c("gene1", "gene1"), 
                   wt1=c(.5, .5), 
                   wt2=c(.25, .25), 
                   stringsAsFactors=FALSE)
  
  expect_warning(prepSNPInfo(si, "Name", "gene"), regexp="Removed 1 duplicate")
  expect_warning(prepSNPInfo(si, "Name", "gene", wt1="wt1"), regexp="Removed 1 duplicate")
  expect_warning(prepSNPInfo(si, "Name", "gene", wt2="wt2"), regexp="Removed 1 duplicate")
  expect_warning(prepSNPInfo(si, "Name", "gene", wt1="wt1", wt2="wt2"), regexp="Removed 1 duplicate")
  expect_warning(prepSNPInfo(si, "Name", "gene", wt1="wt2", wt2="wt1"), regexp="Removed 1 duplicate")
})


test_that("missing snps are dropped", {
  si <- data.frame(Name=c("SNP1", NA), 
                   gene=c("gene1", "gene1"), 
                   wt1=c(.5, .5), 
                   wt2=c(.25, .25), 
                   stringsAsFactors=FALSE)
  
  expect_warning(prepSNPInfo(si, "Name", "gene"), regexp="Removed 1 missing")
  expect_warning(prepSNPInfo(si, "Name", "gene", wt1="wt1"), regexp="Removed 1 missing")
  expect_warning(prepSNPInfo(si, "Name", "gene", wt2="wt2"), regexp="Removed 1 missing")
  expect_warning(prepSNPInfo(si, "Name", "gene", wt1="wt1", wt2="wt2"), regexp="Removed 1 missing")
  expect_warning(prepSNPInfo(si, "Name", "gene", wt1="wt2", wt2="wt1"), regexp="Removed 1 missing")
})


test_that("missing genes are dropped", {
  si <- data.frame(Name=c("SNP1", "SNP2"), 
                   gene=c("gene1", NA), 
                   wt1=c(.5, .5), 
                   wt2=c(.25, .25), 
                   stringsAsFactors=FALSE)
  
  expect_warning(prepSNPInfo(si, "Name", "gene"), regexp="Removed 1 missing")
  expect_warning(prepSNPInfo(si, "Name", "gene", wt1="wt1"), regexp="Removed 1 missing")
  expect_warning(prepSNPInfo(si, "Name", "gene", wt2="wt2"), regexp="Removed 1 missing")
  expect_warning(prepSNPInfo(si, "Name", "gene", wt1="wt1", wt2="wt2"), regexp="Removed 1 missing")
  expect_warning(prepSNPInfo(si, "Name", "gene", wt1="wt2", wt2="wt1"), regexp="Removed 1 missing")
})
