context("protein.info")

# clear out any alterations to the database
suppressMessages(data(thermo))

# test_that somehow affects capture.output so we set up the problem here
protein <- iprotein(c("CSG_METVO", "CSG_METJA"))
suppressMessages(add.obigt())
basis("CHNOS+")
swap.basis("O2", "H2")
pequil <- capture.output(protein.equil(protein, loga.protein=-3))

test_that("protein.equil() reports values consistent with the literature", {
  # the Astar/RT in the paragraph following Eq. 23, p. 6 of DS11
  # (truncated because of rounding)
  expect_true(any(grepl(c("0\\.435.*1\\.36"), pequil)))
  # the log10 activities of the proteins in the left-hand column of the same page
  expect_true(any(grepl(c("-3\\.256.*-2\\.834"), pequil)))
})

# references
# Dick, J. M. and Shock, E. L. (2011) Calculation of the relative chemical stabilities of proteins 
#   as a function of temperature and redox chemistry in a hot spring. 
#   PLoS ONE 6, e22782. http://dx.doi.org/10.1371/journal.pone.0022782
