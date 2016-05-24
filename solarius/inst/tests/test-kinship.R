
context("kinship matrix")

test_that("two custom kinship matrices", {
  data(dat30)
  
  kin2 <- solarKinship2(dat30)
  kin <- kin2 / 2
  
  mod1 <- solarPolygenic(trait1 ~ 1, dat30, kinship = kin)
  mod2 <- solarPolygenic(trait1 ~ 1, dat30, kinship = kin2)
  
  h2r1 <- with(mod1$vcf, Var[varcomp == "h2r"])
  h2r2 <- with(mod2$vcf, Var[varcomp == "h2r"])

  expect_true(h2r1 > h2r2)
})

