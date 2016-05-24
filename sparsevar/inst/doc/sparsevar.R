## ---- eval=FALSE---------------------------------------------------------
#  install.packages("sparsevar", repos = "http://cran.us.r-project.org")

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("devtools", repos = "http://cran.us.r-project.org")
#  devtools::install_github("svazzole/sparsevar")

## ------------------------------------------------------------------------
suppressMessages(library(sparsevar))

## ---- cache = TRUE-------------------------------------------------------
set.seed(1)
sim <- simulateVAR(N = 20, p = 2)

## ---- cache = TRUE-------------------------------------------------------
est <- estimateVAR(sim$data$series, p = 2, options = list(foldsIDs = TRUE))

## ------------------------------------------------------------------------
plotComparisonVAR(sim, est)

## ---- eval=FALSE---------------------------------------------------------
#  results <- estimateVAR(rets)

## ---- eval=FALSE---------------------------------------------------------
#  results <- estimateVAR(rets, p = 3, penalty = "ENET",
#                         options = list(parallel = TRUE, ncores = 5, alpha = 0.95,
#                                        type.measure = "mae", lambda = "lambda.1se"))

## ---- eval=FALSE---------------------------------------------------------
#  sim <- simulateVAR(N = 100, nobs = 250, rho = 0.75, sparsity = 0.05, method = "normal")

