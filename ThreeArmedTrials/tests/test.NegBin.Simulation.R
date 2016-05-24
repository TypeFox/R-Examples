# library(ThreeArmedTrials)
# rm(list=ls())
# 
# rateRef <- 5.1
# ratePla <- 17.4
# 
# Delta <- 0.8
# shape <- 2
# 
# rateExp <- Delta*rateRef + (1-Delta)*ratePla
# 
# #nExp <- nRef <- nPla <- 100
# nExp <- 75
# nRef <- 50
# nPla <- 25
# 
# 
# pRML <- pML <- pSV <- numeric(1000)
# 
# for(i in 1:500){
#   
#   xExp <- rnbinom(n=nExp, mu=rateExp, size=1/shape)
#   xRef <- rnbinom(n=nRef, mu=rateRef, size=1/shape)
#   xPla <- rnbinom(n=nPla, mu=ratePla, size=1/shape)
#   
#   pRML[i] <- taNegbin.test(xExp=xExp, xRef=xRef, xPla=xPla, Delta=Delta, method='RML')$p.value
#   pML[i] <- taNegbin.test(xExp=xExp, xRef=xRef, xPla=xPla, Delta=Delta, method='ML')$p.value
#   pSV[i] <- taNegbin.test(xExp=xExp, xRef=xRef, xPla=xPla, Delta=Delta, method='SampleVariance')$p.value
#   print(i)
# }
# mean(pRML<=0.025)
# mean(pML<=0.025)
# mean(pSV<=0.025)
