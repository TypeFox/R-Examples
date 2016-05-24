context("other bug fixes")

data(seqMetaExample)

###########
# ISSUE N #
###########
## Template


###########
# ISSUE 3 #
###########
## Binomial imputation bug
test_that("imputation inside and outside the package get the same answer)", {
  Z <- Z1
  Z[1, 1] <- NA
  Zi <- Z
  Zi[1, 1] <- mean(Zi[ , 1], na.rm=TRUE)
  
  # binomial case
  so1b<- prepScores(Z=Z, ybin~sex+bmi, family=binomial(), SNPInfo=SNPInfo, data=pheno1)
  so2b<- prepScores(Z=Zi, ybin~sex+bmi, family=binomial(), SNPInfo=SNPInfo, data=pheno1)
  expect_equal(so1b, so2b)
  
  # gaussian case
  so1g<- prepScores(Z=Z, y~sex+bmi, SNPInfo=SNPInfo, data=pheno1)
  so2g<- prepScores(Z=Zi, y~sex+bmi, SNPInfo=SNPInfo, data=pheno1)
  expect_equal(so1g, so2g)
  
  # survival case
  so1c<- prepCox(Z=Z, Surv(time,status)~strata(sex)+bmi, SNPInfo=SNPInfo, data=pheno1)
  so2c<- prepCox(Z=Zi, Surv(time,status)~strata(sex)+bmi, SNPInfo=SNPInfo, data=pheno1)
  expect_equal(so1c, so2c)
  
  # [TBD: Male test cases] - imputation non-trivial as the implemented imputation is not transitive.
})

test_that("Verify the bug does not affect prepScores2)", {
  Z <- Z1
  Z[1, 1] <- NA
  Zi <- Z
  Zi[1, 1] <- mean(Zi[ , 1], na.rm=TRUE)
  
  # binomial case
  so1b<- prepScores2(Z=Z, ybin~sex+bmi, family="binomial", SNPInfo=SNPInfo, data=pheno1)
  so2b<- prepScores2(Z=Zi, ybin~sex+bmi, family="binomial", SNPInfo=SNPInfo, data=pheno1)
  expect_equal(so1b, so2b)
  
  # gaussian case
  so1g<- prepScores2(Z=Z, y~sex+bmi, SNPInfo=SNPInfo, data=pheno1)
  so2g<- prepScores2(Z=Zi, y~sex+bmi, SNPInfo=SNPInfo, data=pheno1)
  expect_equal(so1g, so2g)
  
  # survival case
  so1c<- prepScores2(Z=Z, Surv(time,status)~strata(sex)+bmi, family="cox", SNPInfo=SNPInfo, data=pheno1)
  so2c<- prepScores2(Z=Zi, Surv(time,status)~strata(sex)+bmi, family="cox", SNPInfo=SNPInfo, data=pheno1)
  expect_equal(so1c, so2c)
  
  # [TBD: Male test cases] - imputation non-trivial as the implemented imputation is not transitive.
})


###########
# ISSUE 2 #
###########
# test-issue2.R

###########
# ISSUE 1 #
###########
# test-issue1.R