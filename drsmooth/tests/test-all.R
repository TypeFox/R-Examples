# Run all tests in the inst/tests directory:

library(testthat)
library(drsmooth)

# test scripts do not seem to recognize loaded libraries:
library(car)
library(multcomp)
library(pgirmess)
library(DTK)
library(mgcv)
library(segmented)

test_package("drsmooth")


