## ---- echo = FALSE, message = FALSE--------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")

## ------------------------------------------------------------------------
library(rddtools)
data(house)

## ------------------------------------------------------------------------
house_rdd <- rdd_data(y=house$y, x=house$x, cutpoint=0)

## ----dataPlot------------------------------------------------------------
summary(house_rdd)
plot(house_rdd)

## ----reg_para------------------------------------------------------------
reg_para <- rdd_reg_lm(rdd_object=house_rdd, order=4)
reg_para

plot(reg_para)

## ----RegPlot-------------------------------------------------------------
bw_ik <- rdd_bw_ik(house_rdd)
reg_nonpara <- rdd_reg_np(rdd_object=house_rdd, bw=bw_ik)
print(reg_nonpara)
plot(x=reg_nonpara)

## ----SensiPlot-----------------------------------------------------------
plotSensi(reg_nonpara, from=0.05, to=1, by=0.1)

## ----placeboPlot---------------------------------------------------------
plotPlacebo(reg_nonpara)

## ----DensPlot------------------------------------------------------------
dens_test(reg_nonpara)

## ------------------------------------------------------------------------
set.seed(123)
n_Lee <- nrow(house)
Z <- data.frame(z1 = rnorm(n_Lee, sd=2), 
                z2 = rnorm(n_Lee, mean = ifelse(house<0, 5, 8)), 
                z3 = sample(letters, size = n_Lee, replace = TRUE))
house_rdd_Z <- rdd_data(y = house$y, x = house$x, covar = Z, cutpoint = 0)

## ------------------------------------------------------------------------
## test for equality of means around cutoff:
covarTest_mean(house_rdd_Z, bw=0.3)

## Can also use function covarTest_dis() for Kolmogorov-Smirnov test:
covarTest_dis(house_rdd_Z, bw=0.3)

