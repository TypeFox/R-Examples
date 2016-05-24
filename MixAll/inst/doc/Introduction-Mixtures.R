### R code from vignette source 'Introduction-Mixtures.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: prelim
###################################################
library(MixAll)
MixAll.version <- packageDescription("MixAll")$Version
MixAll.date <- packageDescription("MixAll")$Date


###################################################
### code chunk number 2: Introduction-Mixtures.Rnw:353-355
###################################################
clusterAlgo()
clusterAlgo(algo="SemiSEM",nbIteration=100,epsilon=1e-08)


###################################################
### code chunk number 3: Introduction-Mixtures.Rnw:383-385
###################################################
clusterInit()
clusterInit(method="random", nbInit= 2, algo="CEM", nbIteration=10,epsilon=1e-04)


###################################################
### code chunk number 4: Introduction-Mixtures.Rnw:464-465
###################################################
clusterStrategy()


###################################################
### code chunk number 5: Introduction-Mixtures.Rnw:481-482
###################################################
clusterFastStrategy()


###################################################
### code chunk number 6: Introduction-Mixtures.Rnw:489-490
###################################################
clusterSemiSEMStrategy()


###################################################
### code chunk number 7: Introduction-Mixtures.Rnw:543-546
###################################################
clusterDiagGaussianNames()
clusterDiagGaussianNames("all", "equal", "free")
clusterValidDiagGaussianNames(c("gaussian_pk_sjk","gaussian_p_ljk"))


###################################################
### code chunk number 8: Introduction-Mixtures.Rnw:573-576
###################################################
clusterCategoricalNames()
clusterCategoricalNames("all", "equal")
clusterValidCategoricalNames(c("categorical_pk_pjk","categorical_p_pk"))


###################################################
### code chunk number 9: Introduction-Mixtures.Rnw:613-616
###################################################
clusterPoissonNames()
clusterPoissonNames("all","proportional")
clusterValidPoissonNames(c("poisson_pk_ljk","poisson_p_ljlk"))


###################################################
### code chunk number 10: Introduction-Mixtures.Rnw:689-692
###################################################
clusterGammaNames()
clusterGammaNames("all", "equal","free","free","all")
clusterValidGammaNames(c("gamma_pk_aj_bk","gamma_p_ajk_bjk"))


###################################################
### code chunk number 11: Introduction-Mixtures.Rnw:751-753
###################################################
clusterKernelNames()
clusterValidKernelNames(c("kernelGaussian_pk_sk","kernelGaussian_pk_s"))


###################################################
### code chunk number 12: Introduction-Mixtures.Rnw:861-870
###################################################
data(geyser);
x = as.matrix(geyser); n <- nrow(x); p <- ncol(x);
# add missing values at random
indexes  <- matrix(c(round(runif(10,1,n)), round(runif(10,1,p))), ncol=2);
x[indexes] <- NA;
model <- clusterDiagGaussian(data=x, nbCluster=3, strategy = clusterFastStrategy())
summary(model)
missingValues(model)
plot(model)


###################################################
### code chunk number 13: Introduction-Mixtures.Rnw:904-912
###################################################
data(birds)
x = as.matrix(birds);  n <- nrow(x); p <- ncol(x);
indexes  <- matrix(c(round(runif(10,1,n)), round(runif(10,1,p))), ncol=2);
x[indexes] <- NA;
model <- clusterCategorical(data=x, nbCluster=3, strategy = clusterFastStrategy())
summary(model)
missingValues(model)
plot(model)


###################################################
### code chunk number 14: Introduction-Mixtures.Rnw:944-952
###################################################
data(geyser);
x = as.matrix(geyser); n <- nrow(x); p <- ncol(x);
indexes  <- matrix(c(round(runif(10,1,n)), round(runif(10,1,p))), ncol=2);
x[indexes] <- NA;
model <- clusterGamma(data=x, nbCluster=3, strategy = clusterFastStrategy())
summary(model)
missingValues(model)
plot(model)


###################################################
### code chunk number 15: Introduction-Mixtures.Rnw:983-989
###################################################
data(DebTrivedi)
dt <- DebTrivedi[1:500, c(1, 6,8, 15)]
model <- clusterPoisson( data=dt, nbCluster=3, strategy = clusterFastStrategy())
summary(model)
missingValues(model)
plot(model)


###################################################
### code chunk number 16: Introduction-Mixtures.Rnw:1026-1030
###################################################
data(bullsEye)
model <- clusterKernel( data=bullsEye[,1:2], nbCluster=2, modelNames = "kernelGaussian_pk_s", strategy = clusterFastStrategy())
summary(model)
plot(model)


###################################################
### code chunk number 17: Introduction-Mixtures.Rnw:1059-1067
###################################################
data(HeartDisease.cat)
data(HeartDisease.cont)
ldata = list(HeartDisease.cat,HeartDisease.cont);
lnames = c("categorical_pk_pjk","gaussian_pk_sjk")
model <- clusterHeterogeneous(ldata, lnames, nbCluster=3, strategy = clusterFastStrategy())
summary(model)
missingValues(model)
plot(model)


