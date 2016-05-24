test_that("kinship: kinship without family structure", {
  data(dat30)
  
  kin2 <- solarKinship2(dat30)
  kin <- kin2 / 2

  dat30.subset <- subset(dat30, select = c("id", "trait1", "sex"))
  dat30.subset <- mutate(dat30.subset,
    famid = 0, fa = 0, mo = 0)
  #dat30.subset <- dat30
    
  mod1 <- solarPolygenic(trait1 ~ 1, dat30)
  mod2 <- solarPolygenic(trait1 ~ 1, dat30.subset, kinship = kin)
  
  h2r1 <- with(mod1$vcf, Var[varcomp == "h2r"])
  h2r2 <- with(mod2$vcf, Var[varcomp == "h2r"])

  #expect_equal(h2r1, h2r2)
})

test_that("solarAssoc: different results for different kinship matrices", {
  data(dat50)
  snps <- 1:2
  
  mod1 <- solarAssoc(trait ~ 1, phenodata, snpdata = genodata[, snps], kinship = kin)
  mod2 <- solarAssoc(trait ~ 1, phenodata, snpcovdata = genocovdata[, snps], kinship = 2*kin)
  
  #expect_true(all((mod1$snpf$pSNP != mod2$snpf$pSNP)))
})
