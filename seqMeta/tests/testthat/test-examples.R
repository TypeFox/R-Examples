context("Examples")

data(seqMetaExample)


test_that("basic examples run)", {

  # gaussian case
  cohort1 <- prepScores(Z=Z1, y~sex+bmi, SNPInfo=SNPInfo, data=pheno1)
  cohort2 <- prepScores2(Z=Z1, y~sex+bmi, SNPInfo=SNPInfo, data=pheno1)
  expect_equal(cohort1, cohort2)
  

  cohort1 <- prepScores(Z=Z2, y~sex+bmi, SNPInfo=SNPInfo, kins=kins, data=pheno2)
  cohort2 <- prepScores2(Z=Z2, y~sex+bmi, SNPInfo=SNPInfo, kins=kins, data=pheno2)
  expect_equal(cohort1, cohort2)
  
  
  # binomial case
  cohort1 <- prepScores(Z=Z1, ybin~1, family=binomial(), SNPInfo=SNPInfo, data=pheno1)
  cohort2 <- prepScores2(Z=Z1, ybin~1, family="binomial", SNPInfo=SNPInfo, data=pheno1)  
  expect_equal(cohort1, cohort2)
  
  # survival case
  cohort1 <- prepCox(Z=Z1, Surv(time,status)~strata(sex)+bmi, SNPInfo=SNPInfo, data=pheno1)
  cohort2 <- prepScores2(Z=Z1, Surv(time,status)~strata(sex)+bmi, family="cox", SNPInfo=SNPInfo, data=pheno1)
  expect_equivalent(cohort1, cohort2)
})



test_that("Verify examples run with non-standard names)", {
  si <- SNPInfo
  rm(SNPInfo)
  # gaussian case
  cohort1 <- prepScores(Z=Z1, y~sex+bmi, SNPInfo=si, data=pheno1)
  cohort2 <- prepScores2(Z=Z1, y~sex+bmi, SNPInfo=si, data=pheno1)
  expect_equal(cohort1, cohort2)
  
  
  cohort1 <- prepScores(Z=Z2, y~sex+bmi, SNPInfo=si, kins=kins, data=pheno2)
  cohort2 <- prepScores2(Z=Z2, y~sex+bmi, SNPInfo=si, kins=kins, data=pheno2)
  expect_equal(cohort1, cohort2)
  
  
  # binomial case
  cohort1 <- prepScores(Z=Z1, ybin~1, family=binomial(), SNPInfo=si, data=pheno1)
  cohort2 <- prepScores2(Z=Z1, ybin~1, family="binomial", SNPInfo=si, data=pheno1)  
  expect_equal(cohort1, cohort2)
  
  # survival case
  cohort1 <- prepCox(Z=Z1, Surv(time,status)~strata(sex)+bmi, SNPInfo=si, data=pheno1)
  cohort2 <- prepScores2(Z=Z1, Surv(time,status)~strata(sex)+bmi, family="cox", SNPInfo=si, data=pheno1)
  expect_equivalent(cohort1, cohort2)
})
