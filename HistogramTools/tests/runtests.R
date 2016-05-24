# Simplest possible RUnit boilerplate for now.

library("RUnit")
library("HistogramTools")

# Make the test deterministic - we use a lot of randomly created histograms.
set.seed(0)

# Define tests
testSuite <- defineTestSuite(
  name="HistogramTools Unit Tests",
  dirs=system.file("unitTests", package = "HistogramTools"),
  testFuncRegexp = "^[Tt]est.+")

# Run tests
tests <- runTestSuite(testSuite)

# Print results
printTextProtocol(tests)

# Return success or failure to R CMD CHECK
if (getErrors(tests)$nFail > 0) {
  stop("TEST FAILED!")
}
