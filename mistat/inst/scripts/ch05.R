###################################################
### Chap05Start
###################################################
library(mistat)


###################################################
### ScatterplotChap03
###################################################
data(PLACE)

plot(PLACE[, c("xDev","yDev")],
     main="", xlab="xDev", ylab="yDev")

grid()


###################################################
### MultipleScatterplotChap03
###################################################
plot(PLACE[,-1])


###################################################
### HistogramChap03
###################################################
hist(PLACE[,"tDev"],
     main="", xlab="tDev")


###################################################
### ThreeDimScaterChap03
###################################################
library(lattice)

cloud(tDev ~ xDev*yDev, data=PLACE, type=c("h", "p"), col=1)


###################################################
### BoxPlotChap03
###################################################
boxplot(xDev ~ crcBrd, data=PLACE, 
        ylab="xDev", xlab="crcBrd")


###################################################
### ScatterGroupChap03
###################################################
PLACE$code <- factor(c(rep("lDev", 9*16),
                       rep("mDev", 3*16),
                       rep("hDev", 14*16)))
plot(PLACE[,"xDev"], PLACE[,"yDev"],
     pch=as.integer(PLACE[,"code"]),
     main="", xlab="xDev", ylab="yDev")
grid()
rm(PLACE)


###################################################
### BoxPlotHadpas
###################################################
data(HADPAS)
HADPAS$res7fac <- cut(x=HADPAS$res7, breaks=seq(1300, 2300, by=200), dig.lab=4)
boxplot(res3 ~ res7fac, data = HADPAS)
rm(HADPAS)


###################################################
### MultipleScatterplotALMPIN
###################################################
data(ALMPIN)

plot(ALMPIN)


###################################################
### PlotRelationSocell
###################################################
data(SOCELL)

plot(SOCELL$t1, SOCELL$t2,
     xlab="t1", ylab="t2")


###################################################
### RegressionISC
###################################################
data(SOCELL)

LmISC <- lm(t2 ~ 1 + t1, 
            data=SOCELL)

summary(LmISC)


###################################################
### PlotResidualSocell
###################################################
plot(LmISC, 1)


###################################################
### PlotPredictionIntervalsHadpas
###################################################
data(HADPAS)

plot(HADPAS$res7, HADPAS$res3, xlab="res7", ylab="res3")

LmHadpas <- lm(res3 ~ res7, data=HADPAS)

abline(LmHadpas)

Newdata <- data.frame(res7 =seq(1300, 2300, length.out=200))

Pred <- predict(LmHadpas, 
                newdata=Newdata, 
                interval="predict")

lines(Newdata$res7, Pred[,2])

lines(Newdata$res7, Pred[,3])

rm(Pred, Newdata, LmHadpas)


###################################################
### GasolRegression
###################################################
data(GASOL)

LmYield <- lm(yield ~ 1 + astm + endPt, 
              data=GASOL)

summary(LmYield)


###################################################
### ScatterplotResidualsYield
###################################################
plot(LmYield, 1)


###################################################
### PartialRegressionGasol
###################################################
data(GASOL)

Lm1 <- lm(yield ~ astm, data = GASOL)

Lm2 <- lm(endPt ~ astm, data = GASOL)

plot(Lm2$residuals, Lm1$residuals, 
     xlab=expression(e[2]), 
     ylab=expression(e[1]))


###################################################
### PartialRegression
###################################################
summary(LmYield <- 
          update(object=LmYield, 
                 formula.=. ~ 1 + astm))

summary(LmYield2 <- 
          lm(endPt ~ 1 + astm, 
             data=GASOL))


###################################################
### AlmpinRegression
###################################################
data(ALMPIN)

LmAlmpin <- lm(
  capDiam ~ 1 + diam1 + diam2 + diam3,
  data=ALMPIN)

summary(LmAlmpin)

summary(aov(LmAlmpin))


###################################################
### PlotResidualsAlmpin
###################################################
plot(LmAlmpin, 1, col=1)


###################################################
### GasolStepwiseRegression
###################################################
LmYield <- lm(yield ~ 1 , data=GASOL)

Step <- step(LmYield, direction="both", 
             scope=list(
               lower= ~ 1, 
               upper= ~ endPt + astm + x1 + x2),
             trace=FALSE)

Step$anova


###################################################
### IscRegressionDiagnostic
###################################################
print(influence.measures(LmISC), digits=2)


###################################################
### PlotLeverageLmISC
###################################################
plot(LmISC, 5)


###################################################
### PlotCookDistanceLmAlmpin
###################################################
plot(LmAlmpin, 6)


###################################################
### VendorAnova
###################################################
data(VENDOR)

VENDOR <- stack(VENDOR)               

VENDOR$ind <- as.factor(VENDOR$ind)   

VENDOR$values <- sqrt(VENDOR$values)  

oneway.test(values ~ ind,             
            data = VENDOR, 
            var.equal=T)              

confint(lm(values ~ -1 + ind,         
           data=VENDOR))


###################################################
### PlotBoxplotVendor
###################################################
boxplot(values ~  ind, data=VENDOR)


###################################################
### PlotBoxplotHadpas
###################################################
boxplot(res3 ~ hyb, data=HADPAS)


###################################################
### TukeyHSD
###################################################
HADPAS$hyb <- factor(HADPAS$hyb)

TukeyHSD(aov(res3 ~  hyb, data=HADPAS))


###################################################
### PlotComponentErrorRatesBar
###################################################
data(INSERTION)

barplot(INSERTION$fail / 
          (INSERTION$fail + INSERTION$succ) * 
          100, 
        names.arg=INSERTION$comp, 
        ylab= "%")


###################################################
### Table
###################################################
data(CAR)                     

with(data=CAR,                
     expr=table(cyl, origin)) 


###################################################
### ChisqTest
###################################################
chisq.test(x=CAR$origin, y=CAR$cyl)


###################################################
### Chap05End
###################################################
rm(ALMPIN, CAR, GASOL, HADPAS, INSERTION, SOCELL, VENDOR)
rm(list=ls(pattern="Lm"))
rm(Step)
detach(package:lattice)
detach(package:mistat)
