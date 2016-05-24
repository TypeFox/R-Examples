context("makeup")

test_that("chemical formulas with unknown elements cause a warning", {
  expect_warning(makeup("X"), "element\\(s\\) not in thermo\\$element")
})

test_that("unparseable chemical formulas cause an error", {
  expect_error(makeup("h"), "is not a simple chemical formula")
  expect_error(makeup("1H"), "is not a simple chemical formula")
  expect_error(makeup("(H2O"), "has unpaired parentheses")
  expect_error(makeup("H)2O"), "has unpaired parentheses")
})

test_that("numeric species indices, and coefficients indicating charge can be parsed", {
  # these are all equivalent formulas for the electron
  expect_equal(makeup("-1"), makeup("Z0-1"))
  expect_equal(makeup("-1"), makeup("Z-1+0"))
  # the species index of the electron in thermo$obigt
  ie <- info("e-")
  expect_equal(makeup("-1"), makeup(ie))
})

test_that("signed and fractional coefficients can be parsed", {
  expect_equal(10*makeup("C10.6N1.6P0.1"), makeup("C106N16P1"))
  expect_equal(as.numeric(makeup("C-1.567")), c(1, -1.567))
  expect_equal(as.numeric(makeup("C-1.567+0")), -1.567)
})

test_that("parenthetical and suffixed subformulas can be parsed", {
  expect_equal(makeup("(H)2O"), makeup("H2O"))
})

test_that("summing and multiply formulas works as expected", {
  ispecies <- info(c("B(OH)3", "sepiolite", "lysine:HCl"))
  # the elemental composition all of those, added together
  msaved <- as.array(c(B=1, C=6, Cl=1, H=32, Mg=4, N=2, O=28, Si=6))
  expect_equal(makeup(ispecies, sum=TRUE), msaved)
  # the elemental composition in a 1:2:-1 ratio
  msaved121 <- as.array(c(B=1, C=-6, Cl=-1, H=16, Mg=8, N=-2, O=47, Si=12))
  expect_equal(makeup(ispecies, c(1,2,-1), sum=TRUE), msaved121)
  expect_error(makeup(ispecies, c(1,2), sum=TRUE), "multiplier does not have .* length = number of formulas")
})

test_that("makeup has a fall-through mechanism for matrices and named objects", {
  # this series of tests mimics situations encountered in residue.info() via species.basis()
  # ultimately we want to be able to use species.basis() for species indices, formulas or makeups
  expect_equal(makeup(makeup("CH4")), makeup("CH4"))
  expect_equal(makeup(list(makeup("CH4"))), list(makeup("CH4")))
  basis("CHNOS")
  # we turn the result into a vector using [1, ] so as to drop row names conditionally placed by species.basis
  expect_equal(species.basis("CH4")[1, ], species.basis(info("CH4"))[1, ])
  expect_equal(species.basis(makeup("CH4"))[1, ], species.basis("CH4")[1, ])
  # a matrix should be turned into a list
  protein <- c("LYSC_CHICK", "RNAS1_BOVIN", "CYC_BOVIN", "MYG_PHYCA", "MYG_HORSE")
  pf <- protein.formula(protein)  # a matrix with elements on the columns
  basis(protein)          # yup, a basis set made of proteins, just for fun
  bmat <- basis.elements()  # a matrix with elements on the columns
  expect_equal(as.array(makeup(pf)[[1]]), makeup(as.chemical.formula(pf)[1]))
  expect_equal(makeup(pf), makeup(bmat))
})

test_that("as.chemical.formula moves charge to the end", {
  mkp <- makeup("Z-1HCO3")
  expect_equal(as.chemical.formula(mkp), "HCO3-1")  # i.e. not -1HCO3
})
