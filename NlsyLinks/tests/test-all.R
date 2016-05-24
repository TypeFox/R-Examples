#Modeled after the R6 testing structure: https://github.com/wch/R6/blob/master/tests/testthat.R
library(testthat)
library(NlsyLinks)

testthat::test_check("NlsyLinks")
# test_results <- devtools::test()

# packages <- utils::installed.packages()
# testthatVersion <- packages[packages[, 1]=="testthat", "Version"]
# message("testthat package version: ", testthatVersion)
# 
# if(  testthatVersion >= "0.8" ) {
# #   testthat::test_check("NlsyLinks")
#   system.time(
#     testthat::test_package("NlsyLinks")
#   )
# }
# rm(packages)

#   devtools::test()
