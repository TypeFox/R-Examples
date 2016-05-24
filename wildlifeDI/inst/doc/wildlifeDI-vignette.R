## ------------------------------------------------------------------------
library(wildlifeDI)
data(deer)
deer

## ------------------------------------------------------------------------
deer37 <- deer[1]
deer37
deer38 <- deer[2]
deer38

## ------------------------------------------------------------------------
checkTO(deer37,deer38)

## ----SIFig---------------------------------------------------------------
library(adehabitatHR)
library(rgeos)
#convert ltraj to SpatialPoints - required for kde
pts37 <- SpatialPoints(ld(deer37)[,1:2])
pts38 <- SpatialPoints(ld(deer38)[,1:2])
#compute kernel UD surface - use default method 
#    for obtaining h parameter
kde37 <- kernelUD(pts37)
kde38 <- kernelUD(pts38)
#extract 95% volume contour for HR analysis
hr37 <- getverticeshr(kde37,95)
hr38 <- getverticeshr(kde38,95)
#plot
plot(hr38)
plot(hr37,border="red",add=T)

## ------------------------------------------------------------------------
#Compute SI index
gArea(gIntersection(hr37,hr38))/gArea(gUnion(hr37,hr38))

## ------------------------------------------------------------------------
deers <- GetSimultaneous(deer37,deer38,tc=7.5*60)
deer37.sim <- deers[1]
deer38.sim <- deers[2]
deer37.sim
deer38.sim

## ------------------------------------------------------------------------
Prox(deer37, deer38, tc=7.5*60, dc=50)

## ----label=ProxFig-------------------------------------------------------
prox.df <- Prox(deer37, deer38, tc=7.5*60, dc=50, local=TRUE)
plot(prox.df$date,prox.df$prox,type="l")

## ------------------------------------------------------------------------
spts <- contacts(deer37,deer38,tc=7.5*60, dc=50)
plot(spts)

## ------------------------------------------------------------------------
Ca(deer37, deer38, tc=7.5*60, dc=50)

## ----label=DonFig--------------------------------------------------------
Don(deer37,deer38, tc=7.5*60, dc=50)

## ------------------------------------------------------------------------
Lixn(deer37, deer38, method='spatial', tc=7.5*60, 
     hr1=hr37, hr2=hr38)

## ------------------------------------------------------------------------
Cs(deer37, deer38, tc=7.5*60)

## ------------------------------------------------------------------------
#compute overlap zone
#install.packages('rgeos')
library(rgeos)
oz <- gIntersection(hr37, hr38)

HAI(deer37, deer38, oz, tc=7.5*60, dc=50)

## ------------------------------------------------------------------------
Cr(deer37, deer38, tc=7.5*60)

## ------------------------------------------------------------------------
DI(deer37, deer38, tc=7.5*60)

## ------------------------------------------------------------------------
#obtain the local di analysis data-frame
di.df <- DI(deer37, deer38, tc=7.5*60, local=TRUE)

## ----label=diFig---------------------------------------------------------
#Examine the temporal dynamics of local di
plot(di.df$date, di.df$di,type="l")

## ----label=diFigsmooth---------------------------------------------------
#Smoothed version of local di
di.df$smooth <- 0
#4 fixes/hour x 6 hours on either side of 12 hour centered window
w <- 4*6 
n <- dim(di.df)[1]   #no. of fixes

for (i in (w+1):(n-1-w)){
  di.temp <- di.df$di[(i-w):(i+w)]
  di.df$smooth[i] <- mean(di.temp,na.rm=T)
  }

plot(di.df$date, di.df$smooth,type="l")

## ------------------------------------------------------------------------
IAB(deer37, deer38, dc=50, tc=7.5*60)

## ------------------------------------------------------------------------
df <- IAB(deer37, deer38, dc=50, tc=7.5*60, local=TRUE)
plot(df$date, df$Iab,type='l')

