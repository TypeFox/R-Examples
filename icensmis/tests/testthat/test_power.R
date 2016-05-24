# context("Power test using simulation")
# 
# beta <- 0.5
# HR <- exp(beta)
# sensitivity <- 0.7
# specificity <- 0.98
# testtimes <- 1:4
# blambda <- 0.03
# survivals <- exp(-blambda * testtimes)
# 
# simpower <- function(n, pmiss, pcensor, design, N = 1000) {
#   reject <- sapply(1:n, function(i) {
#     data <- datasim(N, blambda, testtimes, sensitivity, specificity, betas = beta, twogroup = 0.5,
#                     pmiss = pmiss, pcensor = pcensor, design = design)
#     fit <- icmis(ID, testtime, result, data, sensitivity, specificity, ~group)
#     fit$coefficient$p.value < 0.05
#   })
#   mean(reject)
# }
# 
# calpower <- function(pmiss, pcensor, design, N = 1000) {
#   icpower(HR, sensitivity, specificity, survivals, N, pmiss = pmiss, pcensor = pcensor, design = design)$result$power
# }
# 
# test_that("No missing visit, MCAR", {
#   skip_on_cran()
#   expect_less_than(abs(simpower(5000, 0, 0, "MCAR") - calpower(0, 0, "MCAR")), 0.015) 
# })
# 
# test_that("With missing visit, MCAR", {
#   skip_on_cran()
#   expect_less_than(abs(simpower(1000, 0.1, 0, "MCAR") - calpower(0.1, 0, "MCAR")), 0.015) 
#   expect_less_than(abs(simpower(1000, 0.2, 0, "MCAR") - calpower(0.2, 0, "MCAR")), 0.015) 
#   expect_less_than(abs(simpower(1000, 0.3, 0, "MCAR") - calpower(0.3, 0, "MCAR")), 0.015) 
# })
# 
# test_that("No missing visit, NTFP", {
#   skip_on_cran()
#   expect_less_than(abs(simpower(1000, 0, 0, "NTFP", 2000) - calpower(0, 0, "NTFP", 2000)), 0.015) 
# })
# 
# 
# test_that("With missing visit, NTFP", {
#   skip_on_cran()
#   expect_less_than(abs(simpower(3000, 0.1, 0, "NTFP", 2000) - calpower(0.1, 0, "NTFP", 2000)), 0.01) 
#   expect_less_than(abs(simpower(3000, 0.3, 0, "NTFP", 1000) - calpower(0.3, 0, "NTFP", 1000)), 0.01) 
#   expect_less_than(abs(simpower(3000, 0.5, 0, "NTFP", 1000) - calpower(0.5, 0, "NTFP", 1000)), 0.01) 
# })
# 
# test_that("With missing/censoring visit, MCAR", {
#   skip_on_cran()
#   expect_less_than(abs(simpower(2000, 0, 0.1, "MCAR") - calpower(0, 0.1, "MCAR")), 0.015) 
#   expect_less_than(abs(simpower(2000, 0.1, 0.1, "MCAR") - calpower(0.1, 0.1, "MCAR")), 0.015) 
#   expect_less_than(abs(simpower(2000, 0.2, 0.15, "MCAR") - calpower(0.2, 0.15, "MCAR")), 0.015) 
# })
# 
# test_that("With missing/censoring visit, NTFP", {
#   skip_on_cran()
#   expect_less_than(abs(simpower(2000, 0, 0.1, "NTFP", 1500) - calpower(0, 0.1, "NTFP", 1500)), 0.015) 
#   expect_less_than(abs(simpower(2000, 0.1, 0.1, "NTFP", 1500) - calpower(0.1, 0.1, "NTFP", 1500)), 0.015) 
#   expect_less_than(abs(simpower(2000, 0.2, 0.15, "NTFP", 1500) - calpower(0.2, 0.15, "NTFP", 1500)), 0.015) 
# })
