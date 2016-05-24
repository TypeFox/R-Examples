context("Issue 2")
# Monomorphic snps when maf != 0 handled incorrectly
#
# Fix affects prepScores and prepScoresX
#
# Sanity check on prepCox
# TBD
# - prepScores2 

data(seqMetaExample)


###########
# ISSUE 2 #
###########
## Monomophic SNP w/ maf !=0  (CAF == 1  || all hets) handled incorrectly by singlesnpMeta
test_that("is_monomorphic catches all cases)", {
  set.seed(1)
  n <- 100
  x <- cbind(rep(0, n), rep(1, n), rep(2, n), rep(1,n), rep(NA, n), sample(c(0,1,2,NA), n, TRUE))
  colnames(x) <- c("all0", "all1", "all2", "mean1", "allNa", "rand")
  x[1, "mean1"] <- 0
  x[2, "mean1"] <- 2
  
  ans <- c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
  names(ans) <- colnames(x)
  expect_equal(is_monomorphic(x), ans, check.attributes=FALSE)
})

si <- SNPInfo[SNPInfo$gene %in% "gene1", ]
snps_gene1 <- intersect(colnames(Z1), si$Name)
Zgene1 <- Z1[ , snps_gene1]

monos <- c("1000011", "1000009", "1000001")
Zgene1[ , monos[1] ] <- 0
Zgene1[ , monos[2] ] <- 1
Zgene1[ , monos[3] ] <- 2

ans <- rep(FALSE, ncol(Zgene1))
names(ans) <- colnames(Zgene1)
ans[monos] <- TRUE

## Monomophic SNP w/ maf !=0  (CAF == 1  || all hets) handled incorrectly by singlesnpMeta
test_that("Monomophic snps (all 3 cases) handled correctly - prepScores gaussian)", {

  expect_equal(is_monomorphic(Zgene1), ans)
  
  cohort1 <- prepScores(Z=Zgene1, y~sex+bmi, SNPInfo=si, data=pheno1)

  # check maf
  expect_equal(cohort1$gene1$maf[monos[1]], 0, check.attributes=FALSE)
  expect_equal(cohort1$gene1$maf[monos[2]], 0.5, check.attributes=FALSE)
  expect_equal(cohort1$gene1$maf[monos[3]], 1, check.attributes=FALSE)
  
  # check scores
  expect_true(cohort1$gene1$scores[monos[1]] == 0)
  expect_true(cohort1$gene1$scores[monos[2]] == 0)
  expect_true(cohort1$gene1$scores[monos[3]] == 0)

  # check cov
  expect_true(all(cohort1$gene1$cov[ , monos[1]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[1], ] == 0))
  expect_true(all(cohort1$gene1$cov[ , monos[2]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[2], ] == 0))
  expect_true(all(cohort1$gene1$cov[ , monos[3]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[3], ] == 0))
  
  
  out <- singlesnpMeta(cohort1, SNPInfo=si)
  expect_true(out[out$Name == monos[1], "maf"] == 0)
  expect_true(out[out$Name == monos[1], "caf"] == 0)
  expect_true(is.na(out[out$Name == monos[1], "p"]))
  expect_true(is.na(out[out$Name == monos[1], "beta"]))
  expect_true(is.na(out[out$Name == monos[1], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[1], "se"]))
  expect_true(is.infinite(out[out$Name == monos[1], "se.cohort1"])) 
  
  expect_true(out[out$Name == monos[2], "maf"] == 0.5)
  expect_true(out[out$Name == monos[2], "caf"] == 0.5)
  expect_true(is.na(out[out$Name == monos[2], "p"]))
  expect_true(is.na(out[out$Name == monos[2], "beta"]))
  expect_true(is.na(out[out$Name == monos[2], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[2], "se"]))
  expect_true(is.infinite(out[out$Name == monos[2], "se.cohort1"])) 
  
  expect_true(out[out$Name == monos[3], "maf"] == 0)
  expect_true(out[out$Name == monos[3], "caf"] == 1)
  expect_true(is.na(out[out$Name == monos[3], "p"]))
  expect_true(is.na(out[out$Name == monos[3], "beta"]))
  expect_true(is.na(out[out$Name == monos[3], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[3], "se"]))
  expect_true(is.infinite(out[out$Name == monos[3], "se.cohort1"])) 
  
  # test prepScores2 equivalency
  ps2 <- prepScores2(Z=Zgene1, y~sex+bmi, SNPInfo=si, data=pheno1)
  expect_equal(ps2, cohort1)
})

test_that("Monomophic snps (all 3 cases) handled correctly - prepScores binomial)", {
  cohort1 <- prepScores(Z=Zgene1, ybin~1, family=binomial(), SNPInfo=si, data=pheno1)
  
  # check maf
  expect_equal(cohort1$gene1$maf[monos[1]], 0, check.attributes=FALSE)
  expect_equal(cohort1$gene1$maf[monos[2]], 0.5, check.attributes=FALSE)
  expect_equal(cohort1$gene1$maf[monos[3]], 1, check.attributes=FALSE)
  
  # check scores
  expect_true(cohort1$gene1$scores[monos[1]] == 0)
  expect_true(cohort1$gene1$scores[monos[2]] == 0)
  expect_true(cohort1$gene1$scores[monos[3]] == 0)
  
  # check cov
  expect_true(all(cohort1$gene1$cov[ , monos[1]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[1], ] == 0))
  expect_true(all(cohort1$gene1$cov[ , monos[2]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[2], ] == 0))
  expect_true(all(cohort1$gene1$cov[ , monos[3]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[3], ] == 0))
  
  
  out <- singlesnpMeta(cohort1, SNPInfo=si)
  expect_true(out[out$Name == monos[1], "maf"] == 0)
  expect_true(out[out$Name == monos[1], "caf"] == 0)
  expect_true(is.na(out[out$Name == monos[1], "p"]))
  expect_true(is.na(out[out$Name == monos[1], "beta"]))
  expect_true(is.na(out[out$Name == monos[1], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[1], "se"]))
  expect_true(is.infinite(out[out$Name == monos[1], "se.cohort1"])) 
  
  expect_true(out[out$Name == monos[2], "maf"] == 0.5)
  expect_true(out[out$Name == monos[2], "caf"] == 0.5)
  expect_true(is.na(out[out$Name == monos[2], "p"]))
  expect_true(is.na(out[out$Name == monos[2], "beta"]))
  expect_true(is.na(out[out$Name == monos[2], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[2], "se"]))
  expect_true(is.infinite(out[out$Name == monos[2], "se.cohort1"])) 
  
  expect_true(out[out$Name == monos[3], "maf"] == 0)
  expect_true(out[out$Name == monos[3], "caf"] == 1)
  expect_true(is.na(out[out$Name == monos[3], "p"]))
  expect_true(is.na(out[out$Name == monos[3], "beta"]))
  expect_true(is.na(out[out$Name == monos[3], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[3], "se"]))
  expect_true(is.infinite(out[out$Name == monos[3], "se.cohort1"])) 
  
  # test prepScores2 equivalency
  ps2 <- prepScores2(Z=Zgene1, ybin~1, family="binomial", SNPInfo=si, data=pheno1)
  expect_equal(ps2, cohort1)
})

test_that("Monomophic snps (all 3 cases) handled correctly - prepCox)", {
  cohort1 <- prepCox(Z=Zgene1, Surv(time,status)~strata(sex)+bmi, SNPInfo=si, data=pheno1)
  
  # check maf
  expect_equal(cohort1$gene1$maf[monos[1]], 0, check.attributes=FALSE)
  expect_equal(cohort1$gene1$maf[monos[2]], 0.5, check.attributes=FALSE)
  expect_equal(cohort1$gene1$maf[monos[3]], 1, check.attributes=FALSE)
  
  # check scores
  expect_true(cohort1$gene1$scores[monos[1]] == 0)
  expect_true(cohort1$gene1$scores[monos[2]] == 0)
  expect_true(cohort1$gene1$scores[monos[3]] == 0)
  
  # check cov
  expect_true(all(cohort1$gene1$cov[ , monos[1]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[1], ] == 0))
  expect_true(all(cohort1$gene1$cov[ , monos[2]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[2], ] == 0))
  expect_true(all(cohort1$gene1$cov[ , monos[3]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[3], ] == 0))
  
  
  out <- singlesnpMeta(cohort1, SNPInfo=si)
  expect_true(out[out$Name == monos[1], "maf"] == 0)
  expect_true(out[out$Name == monos[1], "caf"] == 0)
  expect_true(is.na(out[out$Name == monos[1], "p"]))
  expect_true(is.na(out[out$Name == monos[1], "beta"]))
  expect_true(is.na(out[out$Name == monos[1], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[1], "se"]))
  expect_true(is.infinite(out[out$Name == monos[1], "se.cohort1"])) 
  
  expect_true(out[out$Name == monos[2], "maf"] == 0.5)
  expect_true(out[out$Name == monos[2], "caf"] == 0.5)
  expect_true(is.na(out[out$Name == monos[2], "p"]))
  expect_true(is.na(out[out$Name == monos[2], "beta"]))
  expect_true(is.na(out[out$Name == monos[2], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[2], "se"]))
  expect_true(is.infinite(out[out$Name == monos[2], "se.cohort1"])) 
  
  expect_true(out[out$Name == monos[3], "maf"] == 0)
  expect_true(out[out$Name == monos[3], "caf"] == 1)
  expect_true(is.na(out[out$Name == monos[3], "p"]))
  expect_true(is.na(out[out$Name == monos[3], "beta"]))
  expect_true(is.na(out[out$Name == monos[3], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[3], "se"]))
  expect_true(is.infinite(out[out$Name == monos[3], "se.cohort1"])) 
  
  #test prepScores2 equivalency
  ps2 <- prepScores2(Z=Zgene1, Surv(time,status)~strata(sex)+bmi, family="cox", SNPInfo=si, data =pheno1)
  expect_equivalent(ps2, cohort1)
})

test_that("Monomophic snps (all 3 cases) handled correctly - prepScores2 gaussian)", {
  cohort1 <- prepScores2(Z=Zgene1, y~sex+bmi, SNPInfo=si, data=pheno1)
  
  # check maf
  expect_equal(cohort1$gene1$maf[monos[1]], 0, check.attributes=FALSE)
  expect_equal(cohort1$gene1$maf[monos[2]], 0.5, check.attributes=FALSE)
  expect_equal(cohort1$gene1$maf[monos[3]], 1, check.attributes=FALSE)
  
  # check scores
  expect_true(cohort1$gene1$scores[monos[1]] == 0)
  expect_true(cohort1$gene1$scores[monos[2]] == 0)
  expect_true(cohort1$gene1$scores[monos[3]] == 0)
  
  # check cov
  expect_true(all(cohort1$gene1$cov[ , monos[1]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[1], ] == 0))
  expect_true(all(cohort1$gene1$cov[ , monos[2]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[2], ] == 0))
  expect_true(all(cohort1$gene1$cov[ , monos[3]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[3], ] == 0))
  
  
  out <- singlesnpMeta(cohort1, SNPInfo=si)
  expect_true(out[out$Name == monos[1], "maf"] == 0)
  expect_true(out[out$Name == monos[1], "caf"] == 0)
  expect_true(is.na(out[out$Name == monos[1], "p"]))
  expect_true(is.na(out[out$Name == monos[1], "beta"]))
  expect_true(is.na(out[out$Name == monos[1], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[1], "se"]))
  expect_true(is.infinite(out[out$Name == monos[1], "se.cohort1"])) 
  
  expect_true(out[out$Name == monos[2], "maf"] == 0.5)
  expect_true(out[out$Name == monos[2], "caf"] == 0.5)
  expect_true(is.na(out[out$Name == monos[2], "p"]))
  expect_true(is.na(out[out$Name == monos[2], "beta"]))
  expect_true(is.na(out[out$Name == monos[2], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[2], "se"]))
  expect_true(is.infinite(out[out$Name == monos[2], "se.cohort1"])) 
  
  expect_true(out[out$Name == monos[3], "maf"] == 0)
  expect_true(out[out$Name == monos[3], "caf"] == 1)
  expect_true(is.na(out[out$Name == monos[3], "p"]))
  expect_true(is.na(out[out$Name == monos[3], "beta"]))
  expect_true(is.na(out[out$Name == monos[3], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[3], "se"]))
  expect_true(is.infinite(out[out$Name == monos[3], "se.cohort1"])) 
})

test_that("Monomophic snps (all 3 cases) handled correctly - prepScores2 binomial)", {
  cohort1 <- prepScores2(Z=Zgene1, ybin~1, family="binomial", SNPInfo=si, data=pheno1)
  
  # check maf
  expect_equal(cohort1$gene1$maf[monos[1]], 0, check.attributes=FALSE)
  expect_equal(cohort1$gene1$maf[monos[2]], 0.5, check.attributes=FALSE)
  expect_equal(cohort1$gene1$maf[monos[3]], 1, check.attributes=FALSE)
  
  # check scores
  expect_true(cohort1$gene1$scores[monos[1]] == 0)
  expect_true(cohort1$gene1$scores[monos[2]] == 0)
  expect_true(cohort1$gene1$scores[monos[3]] == 0)
  
  # check cov
  expect_true(all(cohort1$gene1$cov[ , monos[1]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[1], ] == 0))
  expect_true(all(cohort1$gene1$cov[ , monos[2]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[2], ] == 0))
  expect_true(all(cohort1$gene1$cov[ , monos[3]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[3], ] == 0))
  
  
  out <- singlesnpMeta(cohort1, SNPInfo=si)
  expect_true(out[out$Name == monos[1], "maf"] == 0)
  expect_true(out[out$Name == monos[1], "caf"] == 0)
  expect_true(is.na(out[out$Name == monos[1], "p"]))
  expect_true(is.na(out[out$Name == monos[1], "beta"]))
  expect_true(is.na(out[out$Name == monos[1], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[1], "se"]))
  expect_true(is.infinite(out[out$Name == monos[1], "se.cohort1"])) 
  
  expect_true(out[out$Name == monos[2], "maf"] == 0.5)
  expect_true(out[out$Name == monos[2], "caf"] == 0.5)
  expect_true(is.na(out[out$Name == monos[2], "p"]))
  expect_true(is.na(out[out$Name == monos[2], "beta"]))
  expect_true(is.na(out[out$Name == monos[2], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[2], "se"]))
  expect_true(is.infinite(out[out$Name == monos[2], "se.cohort1"])) 
  
  expect_true(out[out$Name == monos[3], "maf"] == 0)
  expect_true(out[out$Name == monos[3], "caf"] == 1)
  expect_true(is.na(out[out$Name == monos[3], "p"]))
  expect_true(is.na(out[out$Name == monos[3], "beta"]))
  expect_true(is.na(out[out$Name == monos[3], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[3], "se"]))
  expect_true(is.infinite(out[out$Name == monos[3], "se.cohort1"])) 
})

test_that("Monomophic snps (all 3 cases) handled correctly - prepScores2 survival)", {
  cohort1 <- prepScores2(Z=Zgene1, Surv(time,status)~strata(sex)+bmi, family="cox", SNPInfo=si, data =pheno1)
  
  # check maf
  expect_equal(cohort1$gene1$maf[monos[1]], 0, check.attributes=FALSE)
  expect_equal(cohort1$gene1$maf[monos[2]], 0.5, check.attributes=FALSE)
  expect_equal(cohort1$gene1$maf[monos[3]], 1, check.attributes=FALSE)
  
  # check scores
  expect_true(cohort1$gene1$scores[monos[1]] == 0)
  expect_true(cohort1$gene1$scores[monos[2]] == 0)
  expect_true(cohort1$gene1$scores[monos[3]] == 0)
  
  # check cov
  expect_true(all(cohort1$gene1$cov[ , monos[1]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[1], ] == 0))
  expect_true(all(cohort1$gene1$cov[ , monos[2]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[2], ] == 0))
  expect_true(all(cohort1$gene1$cov[ , monos[3]] == 0))
  expect_true(all(cohort1$gene1$cov[ monos[3], ] == 0))
  
  
  out <- singlesnpMeta(cohort1, SNPInfo=si)
  expect_true(out[out$Name == monos[1], "maf"] == 0)
  expect_true(out[out$Name == monos[1], "caf"] == 0)
  expect_true(is.na(out[out$Name == monos[1], "p"]))
  expect_true(is.na(out[out$Name == monos[1], "beta"]))
  expect_true(is.na(out[out$Name == monos[1], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[1], "se"]))
  expect_true(is.infinite(out[out$Name == monos[1], "se.cohort1"])) 
  
  expect_true(out[out$Name == monos[2], "maf"] == 0.5)
  expect_true(out[out$Name == monos[2], "caf"] == 0.5)
  expect_true(is.na(out[out$Name == monos[2], "p"]))
  expect_true(is.na(out[out$Name == monos[2], "beta"]))
  expect_true(is.na(out[out$Name == monos[2], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[2], "se"]))
  expect_true(is.infinite(out[out$Name == monos[2], "se.cohort1"])) 
  
  expect_true(out[out$Name == monos[3], "maf"] == 0)
  expect_true(out[out$Name == monos[3], "caf"] == 1)
  expect_true(is.na(out[out$Name == monos[3], "p"]))
  expect_true(is.na(out[out$Name == monos[3], "beta"]))
  expect_true(is.na(out[out$Name == monos[3], "beta.cohort1"]))
  expect_true(is.na(out[out$Name == monos[3], "se"]))
  expect_true(is.infinite(out[out$Name == monos[3], "se.cohort1"])) 
})

