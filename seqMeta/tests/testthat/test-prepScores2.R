context("prepScores2")

data(seqMetaExample)


test_that("prepScores2 equals prepScores (gaussian w/o family)", {
  ps <- prepScores(Z=Z1, y~sex+bmi, SNPInfo=SNPInfo, data=pheno1)
  ps2 <- prepScores2(Z=Z1, y~sex+bmi, family="gaussian", SNPInfo=SNPInfo, data=pheno1)
  
  expect_equal(ps2, ps)
})


test_that("prepScores2 equals prepScores (gaussian w/ family)", {
  ps <- prepScores(Z=Z2, y~sex+bmi, SNPInfo = SNPInfo, kins=kins, data=pheno2)
  ps2 <- prepScores2(Z=Z2, y~sex+bmi, family="gaussian", SNPInfo=SNPInfo, kins=kins, data=pheno2)
  
  expect_equal(ps2, ps)
})

test_that("prepScores2 equals prepScores (binomial)", {
  ps <- prepScores(Z=Z1, ybin~1, family=binomial(), SNPInfo=SNPInfo, data=pheno1)
  ps2 <- prepScores2(Z=Z1, ybin~1, family="binomial", SNPInfo=SNPInfo, data=pheno1)
  
  expect_equal(ps2, ps)
})

# ###############
# ## SAME TESTS but for X
# ###############
test_that("prepScores2 equals prepScoresX (gaussian w/o family)", {
  ps <- prepScoresX(Z=Z1, y~sex+bmi, male=pheno1$sex-1, SNPInfo=SNPInfo, data=pheno1)
  ps2 <- prepScores2(Z=Z1, y~sex+bmi, male=pheno1$sex-1, family="gaussian", SNPInfo=SNPInfo, data=pheno1)
  
  expect_equal(ps2, ps)
})

test_that("prepScores2 equals prepScoresX (gaussian w/ family)", {
 ps <- prepScoresX(Z=Z2, y~sex+bmi, male=pheno2$sex-1, SNPInfo = SNPInfo, kins=kins, data=pheno2)
 ps2 <- prepScores2(Z=Z2, y~sex+bmi, male=pheno2$sex-1, family="gaussian", SNPInfo=SNPInfo, kins=kins, data=pheno2)
  
 expect_equal(ps2, ps)
})

test_that("prepScores2 equals prepScoresX (binomial)", {
  ps <- prepScoresX(Z=Z1, ybin~1, male=pheno1$sex-1, family=binomial(), SNPInfo=SNPInfo, data=pheno1)
  ps2 <- prepScores2(Z=Z1, ybin~1, male=pheno1$sex-1, family="binomial", SNPInfo=SNPInfo, data=pheno1)
  
  expect_equal(ps2, ps)
})


###############
## SAME TESTS but for Survival
###############
test_that("prepScores2 equals prepCox", {
  ps <- prepCox(Z=Z1, Surv(time,status)~strata(sex)+bmi, SNPInfo=SNPInfo, data=pheno1)
  ps2 <- prepScores2(Z=Z1, Surv(time,status)~strata(sex)+bmi, family="cox", SNPInfo=SNPInfo, data=pheno1)
  
  #prepCox stores the coxfit call in the attributes because it uses "by" instead of 'tapply'
  expect_equal(ps2, ps, check.attributes=FALSE)
})

test_that("prepScores2 fails with incompatiple arguments", {
  expect_error(prepScores2(Z=Z1, Surv(time,status)~strata(sex)+bmi, family="cox", 
                           SNPInfo=SNPInfo, kins=kins, data=pheno1))
  
  expect_error(prepScores2(Z=Z1, Surv(time,status)~strata(sex)+bmi, family="cox", 
                           SNPInfo=SNPInfo, male=pheno1$sex-1, data=pheno1))
})

###############
## Missing values
###############
# test missing all missing genotypes in a snp gives the same as if the snp was removed from the genotype matrix


