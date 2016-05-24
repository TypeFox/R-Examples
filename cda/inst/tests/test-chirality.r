context("Checking that absence of chirality yields no optical activity")

gold <- epsAu(seq(400, 600, by=100))

cl <- cluster_dimer(d=100, 
              dihedral=0*pi/180, alpha1=20*pi/180, alpha2=0,
              a=35, b=12)

circular <- circular_dichroism_spectrum(cl, gold, averaging="GL", Nquad=100)
CD <- subset(circular, type == "CD")
xsec <- subset(circular, type != "CD")

test_that("test that a dimer with plane of symmetry yields no CD", {
  expect_that(max(CD[['value']]) / max(xsec[['value']]) < 1e-10, is_true())
})