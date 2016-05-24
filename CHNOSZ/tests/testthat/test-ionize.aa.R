context("ionize.aa")

test_that("handling of repeated T, P, pH values is correct", {
  expect_message(expect_error(ionize.aa(T=c(25, 25, 100, 100, 100), P=c(100, 1000, 1000, 1000, 1000))),
                              "18 species at 3 values of T and P \\(wet\\)")
  expect_message(ia.25x10 <- ionize.aa(T=rep(25,10), ret.val="pK"), "18 species at 298.15 K and 1 bar \\(wet\\)")
  # we have ten rows of the same values
  expect_identical(ia.25x10[1, ], apply(ia.25x10, 2, unique))
  # also ten rows of the same values for same pH
  ia.7x10 <- ionize.aa(pH=rep(7, 10), ret.val="pK")
  expect_identical(ia.7x10[1, ], apply(ia.7x10, 2, unique))
  # same behaviour with alpha and aavals
  ia.7x10 <- ionize.aa(pH=rep(7, 10), ret.val="alpha")
  expect_identical(ia.7x10[1, ], apply(ia.7x10, 2, unique))
  ia.7x10 <- ionize.aa(pH=rep(7, 10), ret.val="aavals")
  expect_identical(ia.7x10[1, ], apply(ia.7x10, 2, unique))
})

test_that("charge summations for aa compositions are internally consistent", {
  # one pH value should give the same results as two of the same value
  expect_equal(ionize.aa(thermo$protein[1:10,], pH=7), 
    ionize.aa(thermo$protein[1:10,], pH=c(7,7))[1, , drop=FALSE], check.attributes=FALSE)
  # at pH extremes we should be fully ionized
  iplus <- match(c("His", "Lys", "Arg", "chains"), colnames(thermo$protein))
  iminus <- match(c("Cys", "Asp", "Glu", "Tyr", "chains"), colnames(thermo$protein))
  ia520 <- ionize.aa(thermo$protein[1:10,], pH=c(-5, 20))
  expect_equal(ia520[1, ], rowSums(thermo$protein[1:10, iplus]))
  expect_equal(ia520[2, ], -rowSums(thermo$protein[1:10, iminus]))
})

test_that("charge summations for aa compositions are consistent with literature", {
  # comparison with values for LYSC_CHICK digitized from Fig. 10 of Dick et al., 2006
  # charge at 25, 100, 150 degrees and at 25 degrees with cysteine suppressed (oxidized)
  # at pH 4, 6, 8, 10, 12, 14
  Z.LYSC_CHICK.25 <- c(13.3, 8.5, 5.1, -6.2, -14, -20.6)
  Z.LYSC_CHICK.100 <- c(13.3, 7.8, -3.4, -15.2, -20.9, -20.9)
  Z.LYSC_CHICK.150 <- c(13.3, 7.1, -8.0, -20.0, -20.9, -20.9)
  Z.LYSC_CHICK.25_oxid <- c(13.5, 8.7, 7.6, 1.6, -6.2, -12.9)
  aa <- ip2aa(iprotein("LYSC_CHICK"))
  pH <- c(4, 6, 8, 10, 12, 14)
  # the literature values are significantly different at this tolerance (the following is not TRUE)
  # expect_equal(Z.LYSC_CHICK.25, Z.LYSC_CHICK.100, 1e-1)
  expect_equal(ionize.aa(aa, pH=pH, T=25)[, 1], Z.LYSC_CHICK.25, tolerance=1e-1, check.attributes=FALSE)
  expect_equal(ionize.aa(aa, pH=pH, T=100)[, 1], Z.LYSC_CHICK.100, tolerance=1e-1, check.attributes=FALSE)
  expect_equal(ionize.aa(aa, pH=pH, T=150)[, 1], Z.LYSC_CHICK.150, tolerance=1e-1, check.attributes=FALSE)
  expect_equal(ionize.aa(aa, pH=pH, T=25, suppress.Cys=TRUE)[, 1], Z.LYSC_CHICK.25_oxid, tolerance=1e-2, check.attributes=FALSE)
})

test_that("heat capacity and Gibbs energy of ionization are consistent with literature", {
  # heat capacity (kcal mol-1 K-1) of AMYA_PYRFU at 60, 80, 100, 120, 140 degrees
  # for nonionzed protein and ionized protein at pH 6 and 12
  # digitized from Fig. 11 of Dick et al., 2006
  Cp.AMYA_PYRFU.nonion <- c(41.26, 42.22, 42.85, 43.33, 43.74)
  Cp.AMYA_PYRFU.pH6 <- c(37.44, 37.81, 37.78, 37.41, 36.70)
  Cp.AMYA_PYRFU.pH12 <- c(36.19, 36.59, 36.52, 36.11, 35.30)
  aa <- ip2aa(iprotein("AMYA_PYRFU"))
  Cp.ionization.pH6 <- ionize.aa(aa, "Cp", T=c(60, 80, 100, 120, 140), pH=6)
  Cp.ionization.pH12 <- ionize.aa(aa, "Cp", T=c(60, 80, 100, 120, 140), pH=12)
  # the literature values are significantly different at this tolerance (the following is not TRUE)
  # expect_equal(Cp.AMYA_PYRFU.pH6 - Cp.AMYA_PYRFU.nonion, Cp.AMYA_PYRFU.pH12 - Cp.AMYA_PYRFU.nonion, 1e-2)
  expect_equal(Cp.ionization.pH6[,1], (Cp.AMYA_PYRFU.pH6 - Cp.AMYA_PYRFU.nonion)*1000, tolerance=1e-2, check.attributes=FALSE)
  expect_equal(Cp.ionization.pH12[,1], (Cp.AMYA_PYRFU.pH12 - Cp.AMYA_PYRFU.nonion)*1000, tolerance=1e-2, check.attributes=FALSE)
  # Gibbs energy (Mcal mol-1) of AMY_BACSU at pH 0, 2, 4, 6, 8, 10, 12, 14 at 25 and 100 degrees
  # digitized from Fig. 12 of Dick et al., 2006
  G.AMY_BACSU.25 <- c(-24.9, -24.9, -24.7, -24.5, -24.4, -23.9, -23.5, -23.2)
  G.AMY_BACSU.100 <- c(-26.7, -26.7, -26.4, -26.1, -25.7, -25.1, -24.9, -24.9)
  # calculate the Gibbs energies of the nonionized proteins using the same [Met] parameters as in the paper
  add.obigt()
  G.nonionized <- subcrt("AMY_BACSU", T=c(25, 100))$out[[1]]$G
  aa <- ip2aa(iprotein("AMY_BACSU"))
  G.ionization.25 <- ionize.aa(aa, "G", T=25, pH=seq(0, 14, 2))[,1]
  G.ionization.100 <- ionize.aa(aa, "G", T=100, pH=seq(0, 14, 2))[,1]
  expect_equal(G.nonionized[1] + G.ionization.25, G.AMY_BACSU.25 * 1e6, tolerance=1e-3, check.attributes=FALSE)
  expect_equal(G.nonionized[2] + G.ionization.100, G.AMY_BACSU.100 * 1e6, tolerance=1e-3, check.attributes=FALSE)
  # restore thermo$obigt to original state
  data(thermo)
})

test_that("affinity of ionization is consistent with manual calculations", {
  # the equilibrium constant of ionization of [Cys] at 25 and 100 degres C
  logK.Cys <- subcrt(c("[Cys]", "[Cys-]", "H+"), c(-1, 1, 1), T=c(25, 100))$out$logK
  # the affinity (A/2.303RT) of ionization of [Cys] at pH 7
  A.2303RT.Cys.pH7 <- logK.Cys + 7
  # the equilibrium constant of ionization of [His] at 25 degrees C
  logK.His <- subcrt(c("[His]", "H+", "[His+]"), c(-1, -1, 1), T=25)$out$logK
  # the affinity (A/2.303RT) of ionization of [His] at pH 7 and 14
  A.2303RT.His.pH7_14 <- logK.His - c(7, 14)
  # calculate the affinities at pH 7 at 25 and 100 degrees C using ionize()
  A.ionization.pH7 <- ionize.aa(property="A", T=c(25, 100), ret.val="aavals")
  iCys <- match("[Cys-]", colnames(A.ionization.pH7))
  expect_equal(A.2303RT.Cys.pH7, A.ionization.pH7[, iCys], check.attributes=FALSE)
  # calculate the affinities at pH 7 and 14 at 25 degrees C using ionize()
  A.ionization.pH7_14 <- ionize.aa(property="A", pH=c(7, 14), ret.val="aavals")
  iHis <- match("[His+]", colnames(A.ionization.pH7_14))
  expect_equal(A.2303RT.His.pH7_14, A.ionization.pH7_14[, iHis], check.attributes=FALSE)
  # test whether the additive value for a protein is internally consistent
  alpha <- ionize.aa(ret.val="alpha")
  affinity <- ionize.aa(property="A", ret.val="aavals")
  aa <- ip2aa(iprotein("LYSC_CHICK"))
  iionize <- match(c("Cys", "Asp", "Glu", "His", "Lys", "Arg", "Tyr", "chains", "chains"), colnames(aa))
  # the sum of the affinities of ionization of each ionizable group mutiplied by 
  # their degree of formation and by their frequency in the protein
  A.protein <- sum(alpha * affinity * aa[, iionize])
  expect_equal(ionize.aa(aa, property="A")[1, ], A.protein, check.attributes=FALSE)
})

# references

# Dick, J. M., LaRowe, D. E. and Helgeson, H. C. (2006) 
#   Temperature, pressure, and electrochemical constraints on protein speciation: 
#   Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. 
#   Biogeosciences 3, 311--336. http://dx.doi.org/10.5194/bg-3-311-2006
