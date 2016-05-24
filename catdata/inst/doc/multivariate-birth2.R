### R code from vignette source 'multivariate-birth2.Rnw'

###################################################
### code chunk number 1: multivariate-birth2.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: multivariate-birth2.Rnw:19-22 (eval = FALSE)
###################################################
## library(catdata)
## data(birth)
## attach(birth)


###################################################
### code chunk number 3: multivariate-birth2.Rnw:29-36 (eval = FALSE)
###################################################
## intensive <- rep(0,length(Intensive))
## intensive[Intensive>0] <- 1
## Intensive <- intensive
## 
## previous <- Previous
## previous[previous>1] <- 2
## Previous <- previous


###################################################
### code chunk number 4: multivariate-birth2.Rnw:41-42 (eval = FALSE)
###################################################
## library(gee)


###################################################
### code chunk number 5: multivariate-birth2.Rnw:47-54 (eval = FALSE)
###################################################
## library(VGAM)
## Birth <- as.data.frame(na.omit(cbind(Intensive, Cesarean, Sex, Weight, Previous, 
## AgeMother)))
## detach(birth)
## bivarlogit <- vglm(cbind(Intensive , Cesarean) ~ Weight + AgeMother + 
## as.factor(Sex) + as.factor(Previous), binom2.or(zero=NULL), data=Birth)
## summary(bivarlogit)


###################################################
### code chunk number 6: multivariate-birth2.Rnw:59-82 (eval = FALSE)
###################################################
## n <- dim(Birth)[1]
## ID <- rep(1:n,2)
## 
## InterceptInt <- InterceptCes <- rep(1, 2*n)
## InterceptInt[(n+1):(2*n)] <- InterceptCes[1:n] <- 0
## 
## AgeMotherInt <- AgeMotherCes <- rep(Birth$AgeMother,2)
## AgeMotherInt[(n+1):(2*n)] <- AgeMotherCes[1:n] <- 0
## 
## SexInt <- SexCes <- rep(Birth$Sex,2)
## SexInt[SexInt==1] <- SexCes[SexCes==1] <- 0
## SexInt[SexInt==2] <- SexCes[SexCes==2] <- 1
## SexInt[(n+1):(2*n)] <- SexCes[1:n] <- 0
## 
## PrevBase <- rep(Birth$Previous,2)
## PreviousInt1 <- PreviousCes1 <- PreviousInt2 <- PreviousCes2 <- rep(0, 2*n)
## PreviousInt1[PrevBase==1] <- PreviousCes1[PrevBase==1] <- 1
## PreviousInt2[PrevBase>=2] <- PreviousCes2[PrevBase>=2] <- 1
## PreviousInt1[(n+1):(2*n)] <- PreviousInt2[(n+1):(2*n)] <- PreviousCes1[1:n] <- 
##   PreviousCes2[1:n] <- 0
## 
## WeightInt <- WeightCes <- rep(Birth$Weight,2)
## WeightInt[(n+1):(2*n)] <- WeightCes[1:n] <- 0


###################################################
### code chunk number 7: multivariate-birth2.Rnw:87-91 (eval = FALSE)
###################################################
## GeeDat <- as.data.frame(cbind(ID, InterceptInt, InterceptCes, SexInt , SexCes , 
## WeightInt , WeightCes , PreviousInt1 , PreviousInt2, PreviousCes1, 
## PreviousCes2, AgeMotherInt , AgeMotherCes, Response=
## c(Birth$Intensive, Birth$Cesarean)))


###################################################
### code chunk number 8: multivariate-birth2.Rnw:96-102 (eval = FALSE)
###################################################
## gee1 <- gee (Response ~ -1 + InterceptInt + InterceptCes + WeightInt + WeightCes 
##              + AgeMotherInt + AgeMotherCes + SexInt + SexCes +
## PreviousInt1 + PreviousCes1 + PreviousInt2 + PreviousCes2,
## family=binomial(link=logit), id=ID, data=GeeDat)
## 
## summary(gee1)


###################################################
### code chunk number 9: multivariate-birth2.Rnw:107-124 (eval = FALSE)
###################################################
## coefficients(bivarlogit)[1:2]
## coefficients(gee1)[1:2]
## 
## coefficients(bivarlogit)[4:5]
## coefficients(gee1)[3:4]
## 
## coefficients(bivarlogit)[7:8]
## coefficients(gee1)[5:6]
## 
## coefficients(bivarlogit)[10:11]
## coefficients(gee1)[7:8]
## 
## coefficients(bivarlogit)[13:14]
## coefficients(gee1)[9:10]
## 
## coefficients(bivarlogit)[16:17]
## coefficients(gee1)[11:12]


###################################################
### code chunk number 10: multivariate-birth2.Rnw:127-128 (eval = FALSE)
###################################################
## detach(package:VGAM)


