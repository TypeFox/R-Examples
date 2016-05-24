
test_that("dataprep() throws errors when it should", {

  conc <- c(0.0625, 0.125, 0.25, 0.5, 1)
  ntested <- rep(8, 5)
  nalive <- c(1, 4, 4, 7, 8)

  conc.char <- as.character(conc)
  conc.logical <- as.logical(conc)
  ntested.char <- as.character(ntested)
  ntested.logical <- as.logical(ntested)
  nalive.char <- as.character(nalive)
  nalive.logical <- as.logical(nalive)

  conc.dup <- conc[c(1, 1, 3:5)]

  expect_that(dataprep(dose=conc.char, ntot=ntested, nfx=nalive),
    throws_error())
  expect_that(dataprep(dose=conc.logical, ntot=ntested, nfx=nalive),
    throws_error())
  expect_that(dataprep(dose=conc, ntot=ntested.char, nfx=nalive),
    throws_error())
  expect_that(dataprep(dose=conc, ntot=ntested.logical, nfx=nalive),
    throws_error())
  expect_that(dataprep(dose=conc, ntot=ntested, nfx=nalive.char),
    throws_error())
  expect_that(dataprep(dose=conc, ntot=ntested, nfx=nalive.logical),
    throws_error())

  expect_that(dataprep(dose=conc.dup, ntot=ntested, nfx=nalive),
    throws_error())

})
