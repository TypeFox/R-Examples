################################################################################
# TODO LIST
# TODO: Add test for calculateResultType and other functions.

################################################################################
# CHANGE LOG
# 24.12.2014: Second try, update to thestthat 0.8.1
# 05.12.2013: Updated to thestthat 0.8


# Load testthat package.
library(testthat)

# Run all tests.
test_check(package="strvalidator")

#  # Run manually:
# library(strvalidator)
# library(testthat)
# # Set wd to the dev source files, then run:
# test_dir("tests/testthat")