## ----setup, echo=FALSE---------------------------------------------------
  library("knitr")
  opts_chunk$set(fig.width = 7, fig.height = 4)

## ----inst_cran, eval = F-------------------------------------------------
#  install.packages("gmwm")

## ----inst_github, eval = F-----------------------------------------------
#  # Install dependencies
#  install.packages(c("RcppArmadillo","ggplot2","reshape2","devtools"))
#  
#  # Install the package from github
#  devtools::install_github("SMAC-Group/gmwm")

## ----library_include, message=FALSE--------------------------------------
library(gmwm)

## ----gen_process---------------------------------------------------------
# Set seed for reproducibility
set.seed(1)

# Length of the time series
n = 100000

# Simulate AR(1) + WN
xt = gen.gts(AR1(phi=.99, sigma2 = 0.01) + WN(sigma2=1),n)

## ----wv------------------------------------------------------------------
wv = wvar(xt)
plot(wv)

## ----modelTS-------------------------------------------------------------
TS.mod = AR1() + WN()

## ----GMWM----------------------------------------------------------------
model = gmwm.imu(TS.mod, data = xt)

## ----results, echo = F---------------------------------------------------
results = matrix(c(round(model$estimate,4), c(0.99,0.01,1)), 3, 2)
rownames(results) = rownames(model$estimate)
colnames(results) = c("Estim.","True")
results

## ----GMWMplot------------------------------------------------------------
plot(model)

