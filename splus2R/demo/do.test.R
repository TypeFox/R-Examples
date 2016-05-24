cat("This demo shows an example of using do.test\n")
demoLocation <- system.file(package = "splus2R", "examples.t")

##################################################
# First just do.test("examples.t")
# One test should fail.
##################################################
do.test(demoLocation)

##################################################
# Run verbosely: do.test("examples.t", verbose=TRUE)
# This also shows the tests that pass.
do.test(demoLocation, verbose = TRUE)
