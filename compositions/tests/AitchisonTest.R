
require(compositions)

res <- rAitchison(100,theta=c(1,2,3),sigma=ilrvar2clr(diag(c(0.1,2))))
res2<- rnorm.acomp(100,acomp(c(1,2,3)),ilrvar2clr(diag(c(0.1,2))))
plot(res)
plot(res2,add=TRUE,col="red")

dr = dAitchison(res,theta=c(1,2,3),sigma=ilrvar2clr(diag(c(1,2))))
print(dr)



erg<-AitchisonDistributionIntegrals(c(-1,3,-2),ilrvar2clr(-diag(c(1,2))),grid=60)
print(erg)

(myvar<-with(erg, -1/2*ilrvar2clr(solve(clrvar2ilr(beta)))))
(mymean<-with(erg,myvar%*%theta))

with(erg,myvar-clrVar)
with(erg,mymean-clrMean)


AitchisonDistributionIntegrals


res <- rAitchison(100,theta=c(0.5,1,3),sigma=ilrvar2clr(diag(c(0.1,2))))
plot(res)
AitStats = AitchisonDistributionIntegrals(theta=c(0.5,1,3),sigma=ilrvar2clr(diag(c(0.1,2))))
  plot(clrInv(AitStats$clrMean), add=TRUE, col=2, pch=19)
print(AitStats)
