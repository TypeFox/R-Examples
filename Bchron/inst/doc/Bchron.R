## ----eval=FALSE----------------------------------------------------------
#  install.packages('Bchron')

## ------------------------------------------------------------------------
library(Bchron)

## ----fig.align='center',fig.width=6,fig.height=5-------------------------
ages1 = BchronCalibrate(ages=11553,
                        ageSds=230,
                        calCurves='intcal13',
                        ids='Date-1')
summary(ages1)
plot(ages1)

## ----results='hide'------------------------------------------------------
ages2 = BchronCalibrate(ages=c(3445,11553,7456), 
                        ageSds=c(50,230,110), 
                        calCurves=c('intcal13','intcal13','shcal13'))
summary(ages2)
plot(ages2)

## ----fig.align='center',fig.width=6,fig.height=5-------------------------
ages3 = BchronCalibrate(ages=c(3445,11553), 
                        ageSds=c(50,230), 
                        positions=c(100,150), 
                        calCurves=c('intcal13','normal'))
summary(ages3)
plot(ages3,withPositions=TRUE)

## ------------------------------------------------------------------------
data(Glendalough)
print(Glendalough)

## ----results='hide'------------------------------------------------------
GlenOut = Bchronology(ages=Glendalough$ages,
                      ageSds=Glendalough$ageSds, 
                      calCurves=Glendalough$calCurves,
                      positions=Glendalough$position, 
                      positionThicknesses=Glendalough$thickness,
                      ids=Glendalough$id, 
                      predictPositions=seq(0,1500,by=10))

## ----eval=FALSE----------------------------------------------------------
#  help(Bchronology)

## ----eval=FALSE----------------------------------------------------------
#  summary(GlenOut)

## ------------------------------------------------------------------------
summary(GlenOut, type='convergence')
summary(GlenOut, type='outliers')

## ----fig.align='center',fig.width=6,fig.height=5-------------------------
plot(GlenOut,
     main="Glendalough",
     xlab='Age (cal years BP)',
     ylab='Depth (cm)',
     las=1)

## ----results='hide'------------------------------------------------------
predictAges = predict(GlenOut, 
                      newPositions = c(150,725,1500), 
                      newPositionThicknesses=c(5,0,20))

## ------------------------------------------------------------------------
data(TestChronData)
data(TestRSLData)

## ----eval=FALSE----------------------------------------------------------
#  RSLrun = Bchronology(ages=TestChronData$ages,
#                       ageSds=TestChronData$ageSds,
#                       positions=TestChronData$position,
#                       positionThicknesses=TestChronData$thickness,
#                       ids=TestChronData$id,
#                       calCurves=TestChronData$calCurves,
#                       predictPositions=TestRSLData$Depth)
#  RSLrun2 = BchronRSL(RSLrun,
#                      RSLmean=TestRSLData$RSL,
#                      RSLsd=TestRSLData$Sigma,
#                      degree=3)

## ----fig.align='center',fig.width=6,fig.height=5,eval=FALSE--------------
#  summary(RSLrun2)
#  plot(RSLrun2)

## ----results='hide'------------------------------------------------------
data(Sluggan)
SlugDens = BchronDensity(ages=Sluggan$ages,
                         ageSds=Sluggan$ageSds,
                         calCurves=Sluggan$calCurves,
                         numMix=50)

## ----fig.align='center',fig.width=6,fig.height=5-------------------------
plot(SlugDens,xlab='Age (cal years BP)')

## ----eval=FALSE----------------------------------------------------------
#  SlugDensFast = BchronDensityFast(ages=Sluggan$ages,
#                                   ageSds=Sluggan$ageSds,
#                                   calCurves=Sluggan$calCurves)

