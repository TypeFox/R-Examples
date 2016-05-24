context("Issue 1")
# Duplicated SNP in a snpinfo gene pulls from the genotype matrix twice.
#
# This code affects all functions that take SNPInfo as an argument
#
# TBD
# - prepScores2 prepCox equivalent

data(seqMetaExample)

si <- SNPInfo[SNPInfo$gene %in% "gene1", ]
si_dups <- rbind(si, si[1,])

snps_gene1 <- as.character(intersect(colnames(Z1), si$Name))
Zgene1 <- Z1[ , snps_gene1]

test_that("duplicated SNPS in snpinfo gene only get counted once - prepScores gaussian)", {

  cohort1 <- prepScores(Z=Zgene1, y~sex+bmi, SNPInfo=si, data=pheno1)
  cohort2 <- prepScores(Z=Zgene1, y~sex+bmi, SNPInfo=si_dups, data=pheno1)
  expect_equal(length(cohort1), 1)
  expect_equal(length(cohort2), 1)
  expect_equal(cohort1, cohort2)

  expect_equal(length(cohort1$gene1$scores), 15)
  expect_equal(length(cohort2$gene1$scores), 15)
  expect_equal(length(cohort1$gene1$maf), 15)
  expect_equal(length(cohort2$gene1$maf), 15)
  expect_equal(nrow(cohort1$gene1$cov), 15)
  expect_equal(nrow(cohort2$gene1$cov), 15)
  
  # singlesnpMeta 
  single1 <- singlesnpMeta(cohort1, SNPInfo=si, studyBetas=FALSE)
  single2 <- singlesnpMeta(cohort2, SNPInfo=si_dups, studyBetas=FALSE)
  expect_equal(single1, single2)
  expect_equal(nrow(single1), 15)
  expect_equal(nrow(single2), 15)
  
  # burdenMeta
  b1 <- burdenMeta(cohort1, SNPInfo=si)
  b2 <- burdenMeta(cohort2, SNPInfo=si_dups)
  expect_equal(b1, b2)
  
  # skatMeta
  s1 <- skatMeta(cohort1, SNPInfo=si)
  s2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(s1, s2)
  
  # skatOMeta
  sO1 <- skatMeta(cohort1, SNPInfo=si)
  sO2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(sO1, sO2)
  
  # test prepScores2 equivalency
  ps2 <- prepScores2(Z=Zgene1, y~sex+bmi, SNPInfo=si, data=pheno1)
  expect_equal(length(ps2), 1)
  expect_equal(ps2, cohort1)
})

test_that("duplicated SNPS in snpinfo gene only get counted once - prepScores binomial)", {
  cohort1b <- prepScores(Z=Zgene1, ybin~1, family=binomial(), SNPInfo=si, data=pheno1)
  cohort2b <- prepScores(Z=Zgene1, ybin~1, family=binomial(), SNPInfo=si_dups, data=pheno1)
  expect_equal(length(cohort1b), 1)
  expect_equal(length(cohort2b), 1)
  expect_equal(cohort1b, cohort2b)
  
  expect_equal(length(cohort1b$gene1$scores), 15)
  expect_equal(length(cohort2b$gene1$scores), 15)
  expect_equal(length(cohort1b$gene1$maf), 15)
  expect_equal(length(cohort2b$gene1$maf), 15)
  expect_equal(nrow(cohort1b$gene1$cov), 15)
  expect_equal(nrow(cohort2b$gene1$cov), 15)
  
  # test prepScores2 equivalency
  ps2b <- prepScores2(Z=Zgene1, ybin~1, family="binomial", SNPInfo=si, data=pheno1)
  expect_equal(ps2b, cohort1b)
})

test_that("duplicated SNPS in snpinfo gene only get counted once - prepScoresX gaussian)", {
  cohort1 <- prepScoresX(Z=Zgene1, y~sex+bmi, male=pheno1$sex-1, SNPInfo=si, data=pheno1)
  cohort2 <- prepScoresX(Z=Zgene1, y~sex+bmi, male=pheno1$sex-1, SNPInfo=si_dups, data=pheno1)
  expect_equal(length(cohort1), 1)
  expect_equal(length(cohort2), 1)
  expect_equal(cohort1, cohort2)
  
  expect_equal(length(cohort1$gene1$scores), 15)
  expect_equal(length(cohort2$gene1$scores), 15)
  expect_equal(length(cohort1$gene1$maf), 15)
  expect_equal(length(cohort2$gene1$maf), 15)
  expect_equal(nrow(cohort1$gene1$cov), 15)
  expect_equal(nrow(cohort2$gene1$cov), 15)
  
  # singlesnpMeta 
  single1 <- singlesnpMeta(cohort1, SNPInfo=si, studyBetas=FALSE)
  single2 <- singlesnpMeta(cohort2, SNPInfo=si_dups, studyBetas=FALSE)
  expect_equal(single1, single2)
  expect_equal(nrow(single1), 15)
  expect_equal(nrow(single2), 15)
  
  # burdenMeta
  b1 <- burdenMeta(cohort1, SNPInfo=si)
  b2 <- burdenMeta(cohort2, SNPInfo=si_dups)
  expect_equal(b1, b2)
  
  # skatMeta
  s1 <- skatMeta(cohort1, SNPInfo=si)
  s2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(s1, s2)
  
  # skatOMeta
  sO1 <- skatMeta(cohort1, SNPInfo=si)
  sO2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(sO1, sO2)
  
  # test prepScores2 equivalency
  ps2 <- prepScores2(Z=Zgene1, y~sex+bmi, male=pheno1$sex-1, SNPInfo=si, data=pheno1)
  expect_equal(ps2, cohort1)
})

test_that("duplicated SNPS in snpinfo gene only get counted once - prepScoresX binomial)", {
  cohort1b <- prepScoresX(Z=Zgene1, ybin~1, male=pheno1$sex-1, family=binomial(), SNPInfo=si, data=pheno1)
  cohort2b <- prepScoresX(Z=Zgene1, ybin~1, male=pheno1$sex-1, family=binomial(), SNPInfo=si_dups, data=pheno1)
  expect_equal(length(cohort1b), 1)
  expect_equal(length(cohort2b), 1)
  expect_equal(cohort1b, cohort2b)
  
  expect_equal(length(cohort1b$gene1$scores), 15)
  expect_equal(length(cohort2b$gene1$scores), 15)
  expect_equal(length(cohort1b$gene1$maf), 15)
  expect_equal(length(cohort2b$gene1$maf), 15)
  expect_equal(nrow(cohort1b$gene1$cov), 15)
  expect_equal(nrow(cohort2b$gene1$cov), 15)
  
  # test prepScores2 equivalency
  ps2b <- prepScores2(Z=Zgene1, ybin~1, male=pheno1$sex-1, family="binomial", SNPInfo=si, data=pheno1)
  expect_equal(ps2b, cohort1b)
})

test_that("duplicated SNPS in snpinfo gene only get counted once - prepCox)", {
  cohort1 <- prepCox(Z=Zgene1, Surv(time,status)~strata(sex)+bmi, SNPInfo=si, data=pheno1)
  cohort2 <- prepCox(Z=Zgene1, Surv(time,status)~strata(sex)+bmi, SNPInfo=si_dups, data=pheno1)
  expect_equal(length(cohort1), 1)
  expect_equal(length(cohort2), 1)
  expect_equal(cohort1, cohort2)
  
  expect_equal(length(cohort1$gene1$scores), 15)
  expect_equal(length(cohort2$gene1$scores), 15)
  expect_equal(length(cohort1$gene1$maf), 15)
  expect_equal(length(cohort2$gene1$maf), 15)
  expect_equal(nrow(cohort1$gene1$cov), 15)
  expect_equal(nrow(cohort2$gene1$cov), 15)
  
  # singlesnpMeta 
  single1 <- singlesnpMeta(cohort1, SNPInfo=si, studyBetas=FALSE)
  single2 <- singlesnpMeta(cohort2, SNPInfo=si_dups, studyBetas=FALSE)
  expect_equal(single1, single2)
  expect_equal(nrow(single1), 15)
  expect_equal(nrow(single2), 15)
  
  # burdenMeta
  b1 <- burdenMeta(cohort1, SNPInfo=si)
  b2 <- burdenMeta(cohort2, SNPInfo=si_dups)
  expect_equal(b1, b2)
  
  # skatMeta
  s1 <- skatMeta(cohort1, SNPInfo=si)
  s2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(s1, s2)
  
  # skatOMeta
  sO1 <- skatMeta(cohort1, SNPInfo=si)
  sO2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(sO1, sO2)
  
  # test prepScores2 equivalency
   ps2 <- prepScores2(Z=Zgene1, Surv(time,status)~strata(sex)+bmi, family="cox", SNPInfo=si, data=pheno1)
   expect_equivalent(ps2, cohort1)
})

test_that("duplicated SNPS in snpinfo gene only get counted once - prepScores2 gaussian)", {
  cohort1 <- prepScores2(Z=Zgene1, y~sex+bmi, SNPInfo=si, data=pheno1)
  cohort2 <- prepScores2(Z=Zgene1, y~sex+bmi, SNPInfo=si_dups, data=pheno1)
  expect_equal(length(cohort1), 1)
  expect_equal(length(cohort2), 1)
  expect_equal(cohort1, cohort2)
  
  expect_equal(length(cohort1$gene1$scores), 15)
  expect_equal(length(cohort2$gene1$scores), 15)
  expect_equal(length(cohort1$gene1$maf), 15)
  expect_equal(length(cohort2$gene1$maf), 15)
  expect_equal(nrow(cohort1$gene1$cov), 15)
  expect_equal(nrow(cohort2$gene1$cov), 15)
  
  # singlesnpMeta 
  single1 <- singlesnpMeta(cohort1, SNPInfo=si, studyBetas=FALSE)
  single2 <- singlesnpMeta(cohort2, SNPInfo=si_dups, studyBetas=FALSE)
  expect_equal(single1, single2)
  expect_equal(nrow(single1), 15)
  expect_equal(nrow(single2), 15)
  
  # burdenMeta
  b1 <- burdenMeta(cohort1, SNPInfo=si)
  b2 <- burdenMeta(cohort2, SNPInfo=si_dups)
  expect_equal(b1, b2)
  
  # skatMeta
  s1 <- skatMeta(cohort1, SNPInfo=si)
  s2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(s1, s2)
  
  # skatOMeta
  sO1 <- skatMeta(cohort1, SNPInfo=si)
  sO2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(sO1, sO2)
  
  # test prepScores equivalency
  ps <- prepScores(Z=Zgene1, y~sex+bmi, SNPInfo=si_dups, data=pheno1)
  expect_equal(ps, cohort2)
})
 
test_that("duplicated SNPS in snpinfo gene only get counted once - prepScores2 binomial)", { 
  cohort1b <- prepScores2(Z=Zgene1, ybin~1, family="binomial", SNPInfo=si, data=pheno1)
  cohort2b <- prepScores2(Z=Zgene1, ybin~1, family="binomial", SNPInfo=si_dups, data=pheno1)
  expect_equal(length(cohort1b), 1)
  expect_equal(length(cohort2b), 1)
  expect_equal(cohort1b, cohort2b)
  
  expect_equal(length(cohort1b$gene1$scores), 15)
  expect_equal(length(cohort2b$gene1$scores), 15)
  expect_equal(length(cohort1b$gene1$maf), 15)
  expect_equal(length(cohort2b$gene1$maf), 15)
  expect_equal(nrow(cohort1b$gene1$cov), 15)
  expect_equal(nrow(cohort2b$gene1$cov), 15)
  
  # test prepScores equivalency
  ps <- prepScores(Z=Zgene1, ybin~1, family=binomial(), SNPInfo=si_dups, data=pheno1)
  expect_equal(ps, cohort2b)
})

test_that("duplicated SNPS in snpinfo gene only get counted once - prepScores2 gaussian w/ kinship)", {
  cohort1 <- prepScores2(Z=Zgene1, y~sex+bmi, male=pheno1$sex-1, SNPInfo=si, data=pheno1)
  cohort2 <- prepScores2(Z=Zgene1, y~sex+bmi, male=pheno1$sex-1, SNPInfo=si_dups, data=pheno1)
  expect_equal(cohort1, cohort2)
  expect_equal(length(cohort1$gene1$scores), 15)
  expect_equal(length(cohort2$gene1$scores), 15)
  expect_equal(length(cohort1$gene1$maf), 15)
  expect_equal(length(cohort2$gene1$maf), 15)
  expect_equal(nrow(cohort1$gene1$cov), 15)
  expect_equal(nrow(cohort2$gene1$cov), 15)
  
  # singlesnpMeta 
  single1 <- singlesnpMeta(cohort1, SNPInfo=si, studyBetas=FALSE)
  single2 <- singlesnpMeta(cohort2, SNPInfo=si_dups, studyBetas=FALSE)
  expect_equal(single1, single2)
  expect_equal(nrow(single1), 15)
  expect_equal(nrow(single2), 15)
  
  # burdenMeta
  b1 <- burdenMeta(cohort1, SNPInfo=si)
  b2 <- burdenMeta(cohort2, SNPInfo=si_dups)
  expect_equal(b1, b2)
  
  # skatMeta
  s1 <- skatMeta(cohort1, SNPInfo=si)
  s2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(s1, s2)
  
  # skatOMeta
  sO1 <- skatMeta(cohort1, SNPInfo=si)
  sO2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(sO1, sO2)
  
  # test prepScoresX equivalency
  ps <- prepScoresX(Z=Zgene1, y~sex+bmi, male=pheno1$sex-1, SNPInfo=si_dups, data=pheno1)
  expect_equal(ps, cohort2)
})

test_that("duplicated SNPS in snpinfo gene only get counted once - prepScores2 binomial w/ kinship)", {
  cohort1b <- prepScores2(Z=Zgene1, ybin~1, male=pheno1$sex-1, family="binomial", SNPInfo=si, data=pheno1)
  cohort2b <- prepScores2(Z=Zgene1, ybin~1, male=pheno1$sex-1, family="binomial", SNPInfo=si_dups, data=pheno1)
  expect_equal(length(cohort1b), 1)
  expect_equal(length(cohort2b), 1)
  expect_equal(cohort1b, cohort2b)
  expect_equal(length(cohort1b$gene1$scores), 15)
  expect_equal(length(cohort2b$gene1$scores), 15)
  expect_equal(length(cohort1b$gene1$maf), 15)
  expect_equal(length(cohort2b$gene1$maf), 15)
  expect_equal(nrow(cohort1b$gene1$cov), 15)
  expect_equal(nrow(cohort2b$gene1$cov), 15)
  
  # test prepScoresX equivalency
  ps <- prepScoresX(Z=Zgene1, ybin~1, male=pheno1$sex-1, family=binomial(), SNPInfo=si_dups, data=pheno1)
  expect_equal(ps, cohort2b)
})

test_that("duplicated SNPS in snpinfo gene only get counted once - prepScores2 survival)", {
  cohort1 <- prepScores2(Z=Zgene1, Surv(time,status)~strata(sex)+bmi, family="cox", SNPInfo=si, data=pheno1)
  cohort2 <- prepScores2(Z=Zgene1, Surv(time,status)~strata(sex)+bmi, family="cox", SNPInfo=si_dups, data=pheno1)
  expect_equal(length(cohort1), 1)
  expect_equal(length(cohort2), 1)
  expect_equal(cohort1, cohort2)
  
  expect_equal(length(cohort1$gene1$scores), 15)
  expect_equal(length(cohort2$gene1$scores), 15)
  expect_equal(length(cohort1$gene1$maf), 15)
  expect_equal(length(cohort2$gene1$maf), 15)
  expect_equal(nrow(cohort1$gene1$cov), 15)
  expect_equal(nrow(cohort2$gene1$cov), 15)
  
  # singlesnpMeta 
  single1 <- singlesnpMeta(cohort1, SNPInfo=si, studyBetas=FALSE)
  single2 <- singlesnpMeta(cohort2, SNPInfo=si_dups, studyBetas=FALSE)
  expect_equal(single1, single2)
  expect_equal(nrow(single1), 15)
  expect_equal(nrow(single2), 15)
  
  # burdenMeta
  b1 <- burdenMeta(cohort1, SNPInfo=si)
  b2 <- burdenMeta(cohort2, SNPInfo=si_dups)
  expect_equal(b1, b2)
  
  # skatMeta
  s1 <- skatMeta(cohort1, SNPInfo=si)
  s2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(s1, s2)
  
  # skatOMeta
  sO1 <- skatMeta(cohort1, SNPInfo=si)
  sO2 <- skatMeta(cohort2, SNPInfo=si_dups)
  expect_equal(sO1, sO2)
  
  # test prepScores2 equivalency
  pcox <- prepCox(Z=Zgene1, Surv(time,status)~strata(sex)+bmi, SNPInfo=si, data=pheno1)
  expect_equivalent(pcox, cohort2)
})
