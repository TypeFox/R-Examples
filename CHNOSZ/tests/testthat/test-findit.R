context("findit")

# this is a long test ... skip it if we're on R CMD check --as-cran
if(!any(grepl("R_CHECK_TIMINGS", names(Sys.getenv())))) {

test_that("findit() returns known values encoded in a species distribution", {
  # test the findit function in three dimensions
  # first calculate equilibrium activities of some proteins
  # using known activities of basis species
  ip <- 1:20
  basis("CHNOS")
  basis("CO2", -pi)        # -3.141593
  basis("H2O", -exp(1))    # -2.718282
  basis("NH3", -sqrt(2))   # -1.414214
  a <- affinity(iprotein=ip)
  # scale relative abundances such that total activity of residues 
  # is unity (since loga.balance=0 is default for findit)
  e <- equilibrate(a, loga.balance=0)
  # loga2 are the equilibrium logarithms of activities of the proteins
  loga2 <- as.numeric(e$loga.equil)
  # return to default values for activities of basis species
  basis("CHNOS")
  ## we have diverged from the reference activities of proteins
  #a <- affinity(iprotein=ip)
  #d <- diagram(a, plot.it=FALSE, loga.balance=0)
  #r0 <- revisit(d, "rmsd", loga.ref, main="")
  #title(main=paste("log activities of 20 proteins for basis activities",
  #  "from numerical constants (ref) and CHNOSZ default (calc)",sep="\n"))
  # now find the activities of the basis species
  # that get us close to reference activities of proteins
  f <- findit(lims=list(CO2=c(-5,0), H2O=c(-5,0), NH3=c(-5,0)),
    objective="RMSD", niter=2, iprotein=ip, loga2=loga2, res=24, rat=0.2, plot.it=FALSE)
  # sanity check: the output values are all the same length
  expect_equal(length(unique(sapply(f$value, length))), 1)
  # -pi, -e and -sqrt(2) were approximately retrieved!
  expect_equal(tail(f$value[[1]],1), -pi, tolerance=1e-2)
  expect_equal(tail(f$value[[2]],1), -exp(1), tolerance=1e-2)
  expect_equal(tail(f$value[[3]],1), -sqrt(2), tolerance=1e-1)
  # we could decrease the tolerance by increasing the resolution and/or iterations in findit()
})

}
