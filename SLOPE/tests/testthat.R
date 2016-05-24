library(testthat)
library(SLOPE)

# Set a random seed for reproducibility when running R CMD CHECK.
# Note that this file is not run by devtools::test, so the tests will be truly
# random when run from RStudio.
set.seed(56789)

# Run the test suite.
test_check("SLOPE")

# Reset the random seed.
rm(.Random.seed)
