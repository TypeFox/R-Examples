context("ionize.aa_pK")

test_that("pK values as a function of temperature are consistent with literature values", {
  # pK at 0, 100, 200, 300 degrees C digitized from Fig. 4 of Dick et al., 2006
  DLH06.pK <- list(`[Cys-]` = c(8.82, 7.44, 7.14, 7.90),
                   `[Asp-]` = c(4.06, 3.80, 4.17, 5.30), 
                   `[Glu-]` = c(4.36, 4.39, 5.11, 6.68), 
                   `[His+]` = c(6.44, 4.91, 3.92, 3.17),
                   `[Lys+]` = c(10.92,8.02, 6.27, 4.75),
                   `[Arg+]` = c(13.38,10.02,8.11, 6.56),
                   `[Tyr-]` = c(9.68, 8.02, 7.59, 8.33),
                   `[AABB+]`= c(2.34, 2.36, 2.62, 2.81),
                   `[AABB-]`= c(10.16,7.85, 6.51, 6.28))
  this.pK <- ionize.aa(T=c(0, 100, 200, 300), ret.val="pK")
  # ionization of [Cys] and [His] is off more than the others
  expect_equal(this.pK[, 1], DLH06.pK$`[Cys-]`, tolerance=1e-1, check.attributes=FALSE)
  expect_equal(this.pK[, 2], DLH06.pK$`[Asp-]`, tolerance=1e-2, check.attributes=FALSE)
  expect_equal(this.pK[, 3], DLH06.pK$`[Glu-]`, tolerance=1e-2, check.attributes=FALSE)
  expect_equal(this.pK[, 4], DLH06.pK$`[His+]`, tolerance=1e-1, check.attributes=FALSE)
  expect_equal(this.pK[, 5], DLH06.pK$`[Lys+]`, tolerance=1e-2, check.attributes=FALSE)
  expect_equal(this.pK[, 6], DLH06.pK$`[Arg+]`, tolerance=1e-2, check.attributes=FALSE)
  expect_equal(this.pK[, 7], DLH06.pK$`[Tyr-]`, tolerance=1e-2, check.attributes=FALSE)
  expect_equal(this.pK[, 8], DLH06.pK$`[AABB+]`,tolerance=1e-2, check.attributes=FALSE)
  expect_equal(this.pK[, 9], DLH06.pK$`[AABB-]`,tolerance=1e-2, check.attributes=FALSE)
})

test_that("there is one pK value for each ionizable group at each temperature", {
  pH <- seq(0, 14, 2)
  T <- seq(0, 150, 15)
  val <- expand.grid(pH=pH, T=T)
  Z <- ionize.aa(pH=val$pH, T=val$T, ret.val="pK")
  expect_equal(length(unique(Z)), 99)
})

# references

# Dick, J. M., LaRowe, D. E. and Helgeson, H. C. (2006) 
#   Temperature, pressure, and electrochemical constraints on protein speciation: 
#   Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. 
#   Biogeosciences 3, 311--336. http://dx.doi.org/10.5194/bg-3-311-2006
