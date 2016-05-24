### R code from vignette source 'vignette_HiLMM.Rnw'

###################################################
### code chunk number 1: loadparameters
###################################################
library(HiLMM)


###################################################
### code chunk number 2: HiLMM_R
###################################################
data_sim= data_simu(500,1000,0.6,0.5)
res_herit=estim_herit(data_sim$Y,data_sim$W)


###################################################
### code chunk number 3: resultHiLMM
###################################################
data_sim$W[1:10,1:10]
data_sim$Y[1:10]


###################################################
### code chunk number 4: resultHiLMM
###################################################
res_herit$heritability
res_herit$CI_low
res_herit$CI_up
res_herit$standard_deviation


###################################################
### code chunk number 5: sessionInfo
###################################################
sessionInfo()


