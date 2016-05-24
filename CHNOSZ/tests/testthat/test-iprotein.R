context("iprotein")

# clear out any prior database alterations
suppressMessages(data(thermo))

test_that("basic searches and conversions work as expected", {
  expect_equal(iprotein(c("LYSC_CHICK", "MYGPHYCA")), c(6, NA))
  # factors causing problems again ...
  f <- system.file("extdata/protein/DS11.csv", package="CHNOSZ")
  aa <- read.csv(f)
  # this adds the proteins
  ip <- add.protein(aa)
  # the replaces the proteins (with the same ones)
  expect_error(ip <- add.protein(aa), "converting factors causes problems replacing protein data")
  # ... should use read.csv(file, stringsAsFactors=FALSE)
})

test_that("errors and messages occur in some circumstances", {
  expect_message(iprotein(c("LYSC_CHICK", "MYGPHYCA")), "1 protein not matched")
  expect_error(seq2aa("LYS_CHICK", "XXX"), "no characters match an amino acid")
  expect_error(add.protein(count.aa("AAA")), "not a data frame with the same columns as thermo\\$protein")
  expect_message(add.protein(ip2aa(iprotein("CYC_BOVIN"))), "replaced 1 existing protein\\(s\\)")
})

test_that("group additivity for proteins gives expected values", {
  # values for chicken lysozyme calculated using group additivity values
  # from Dick et al., 2006 (Biogeosciences 3, 311-336)
  G <- -4206050
  Cp <- 6415.5
  V <- 10421
  formula <- "C613H959N193O185S10"
  # use add.obigt to load the parameters for [Met] sidechain group
  # from above reference instead of the updated values
  # from LaRowe and Dick, 2012 (Geochim Cosmochim Acta 80, 70-91)
  add.obigt()
  lprop <- info(info("LYSC_CHICK"))
  expect_equal(G, lprop$G)
  expect_equal(Cp, lprop$Cp, tolerance=1e-5)
  expect_equal(V, lprop$V, tolerance=1e-4)
  expect_equal(formula, lprop$formula)
})

test_that("amino acid counts taken from a fasta file can be added",{
  ffile <- system.file("extdata/fasta/EF-Tu.aln", package="CHNOSZ")
  aa <- read.fasta(ffile)
  expect_message(ip1 <- add.protein(aa), "added 8 new protein\\(s\\)")
  expect_message(ip2 <- add.protein(aa), "replaced 8 existing protein\\(s\\)")
  # add.protein should return the correct indices for existing proteins
  expect_equal(ip1, ip2)
})

# for the future... make info() faster!
# (especially the loading of ionizable groups for proteins)
#test_that("calculations for ionized proteins meet performance expectations", {
#  expect_that({
#    basis("CHNOS+")
#    i <- info(c("LYSC_CHICK","RNAS1_BOVIN"))
#    species(i)
#    a <- affinity()
#  }, takes_less_than(0.4))
#})
