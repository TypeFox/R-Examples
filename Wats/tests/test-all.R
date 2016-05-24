# library(testthat)
library(Wats)

packages <- utils::installed.packages()
testthatVersion <- packages[packages[, 1]=="testthat", "Version"]
message("testthat package version: ", testthatVersion)

if(  testthatVersion >= "0.8" ) {
  testthat::test_check("Wats")
}
rm(packages)
# testthat::test_package("Wats")
