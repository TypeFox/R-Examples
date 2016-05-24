# reproducing the calculations and graphics from:
# Wijngaard, J. B., Klein Tank, A. M. G. and Konnen, G. P. (2003)
# Homogeneity of 20th century European daily temperature and precipitation series. 
# Int. J. Climatol., 23: 679-692. doi: 10.1002/joc.906

setwd(system.file('extdata/', package='climtrends'))

# data from the ECA station series of Eelde (The Netherlands) from 1900 to 2000
EtmgegData<-ReadEtmgegFile('etmgeg_280.txt',c(2,12,13,15))
until2000<-EtmgegData[which(EtmgegData[,1]<'2001-01-01'),] # zr are the data until 2000
until2000[,2:4]<-until2000[,2:4]/10 # divide by 10 to get a scale of 1 degrees Celsius
u2000AnnualMean<-YearMeanFromDay(until2000,1,2) # calculate yearly mean
u2000AnnualMax<-YearMeanFromDay(until2000,1,4) # calculate max of yearly mean
u2000AnnualMin<-YearMeanFromDay(until2000,1,3) # calculate min of yearly mean
# calculate DTR
u2000DTR<-until2000[,-2] #get rid of the mean
u2000DTR[,2]<-u2000DTR[,3]-u2000DTR[,2] #DTR=MAX-MIN
u2000DTR<-u2000DTR[,-3]
u2000mDTR<-YearMeanFromDay(u2000DTR,1,2)
# calculate VDTR
u2000VDTR<-VDTR(u2000DTR,1,2)

# plot Figure 1
original.parameters<-par()
par(mar=c(1,1,1,1),oma=c(2, 4, 1,4))
plot(1907:2000,u2000mDTR[2:95,2], type = "s",xlim=c(1906,2000),ylim=c(5,13), axes = 1:2, xlab = "", ylab = "",main='HOMOGENEITY OF EUROPEAN CLIMATE SERIES') # for stair steps
mtext("mDTR (deg. Celsius)", side=2, line=3)
par(new=T)
plot(1907:2000, u2000AnnualMean[2:95,2], type = "s",xlim=c(1906,2000),ylim=c(3,11),col="grey",axes = FALSE, bty = "n",xlab = "", ylab = "") # for stair steps
mtext("Mean temperature (deg. Celsius)", side=4,col="grey",line=3)
axis(side=4, col="grey",col.axis="grey")
abline(h = 4:10, col = "gray", lty=3)
abline(v = c(1948,1951,1959), lty=2)
par(new=T)
x<-2:95 #  from 1900 to 2000
y<-u2000mDTR[x,2]
y.loess <- loess(y ~ x, span=0.29, data.frame(x=x, y=y))
y.predict <- predict(y.loess, data.frame(x=x))
plot(u2000mDTR[x,1],y.predict, type = "l",xlim=c(1906,2000),ylim=c(5,13),axes = FALSE, bty = "n",xlab = "", ylab = "")
par(new=T)
x<-2:95 #  from 1900 to 2000
y<-u2000AnnualMean[x,2]
y.loess <- loess(y ~ x, span=0.29, data.frame(x=x, y=y))
y.predict <- predict(y.loess, data.frame(x=x))
plot(u2000AnnualMean[x,1],y.predict, type = "l",xlim=c(1906,2000),ylim=c(3,11),col="grey",axes = FALSE, bty = "n",xlab = "", ylab = "")
par(original.parameters)

dev.new()

#fig. 2
data(Buishand.Critical.Values)
data(Pettitt.Critical.Values)
data(SNHT.Critical.Values)
wu<-ts(u2000mDTR[2:95,2],  start=c(u2000mDTR[1]), frequency=1)
SNHTmDTR<-SNHTabsolute(wu)
BuishandRangemDTR<-BuishandRangeTest(wu)
PettittmDTR<-PettittTest(wu)
par(mfrow=c(1,3))
plot(u2000mDTR[2:95,1],c(SNHTmDTR,NA), type = "l",ylab='Value statistic',xlab='',main='SNHT, mDTR')
c1 <- SNHT.Critical.Values["100","99%"]
c5 <- SNHT.Critical.Values["100","95%"]
negC <- NegCurve(SNHTmDTR)
abline(h=negC*c1,lty="dashed")
abline(h=negC*c5,lty="dotted")
plot(u2000mDTR[2:95,1],BuishandRangemDTR, type = "l",ylab='Value statistic',xlab='',main='Buishand range, mDTR')
c1 <- Buishand.Critical.Values["100","R99%"]
c5 <- Buishand.Critical.Values["100","R95%"]
negC <- NegCurve(BuishandRangemDTR)
abline(h=negC*c1,lty="dashed")
abline(h=negC*c5,lty="dotted")
plot(u2000mDTR[2:95,1],PettittmDTR, type = "l",ylab='Value statistic',xlab='',main='Pettitt test, mDTR')
c1 <- Pettitt.Critical.Values["100","99%"]
c5 <- Pettitt.Critical.Values["100","95%"]
negC <- NegCurve(PettittmDTR)
abline(h=negC*c1,lty="dashed")
abline(h=negC*c5,lty="dotted")


dev.new()

# Figure 3
#SNHT applied to the annual mean (left), maximum (middle) and minimum (right) temperature series
par(mfrow=c(1,3))
wu<-ts(u2000AnnualMean[,2],  start=c(min(u2000mDTR[1]), 1), frequency=1)
tSNHT<-SNHTabsolute(wu)
plot(tSNHT, type = "l",main='SNHT, mean temperature', ylim=c(0,30))
c1 <- SNHT.Critical.Values["100","99%"]
c5 <- SNHT.Critical.Values["100","95%"]
negC <- NegCurve(tSNHT)
abline(h=negC*c1,lty="dashed")
abline(h=negC*c5,lty="dotted")
wu<-ts(u2000AnnualMax[,2],  start=c(min(u2000mDTR[1]), 1), frequency=1)
tSNHT<-SNHTabsolute(wu)
plot(tSNHT, type = "l",main='SNHT, max temperature', ylim=c(0,30))
c1 <- SNHT.Critical.Values["100","99%"]
c5 <- SNHT.Critical.Values["100","95%"]
negC <- NegCurve(tSNHT)
abline(h=negC*c1,lty="dashed")
abline(h=negC*c5,lty="dotted")
wu<-ts(u2000AnnualMin[,2],  start=c(min(u2000mDTR[1]), 1), frequency=1)
tSNHT<-SNHTabsolute(wu)
plot(tSNHT, type = "l",main='SNHT, min temperature', ylim=c(0,30))
c1 <- SNHT.Critical.Values["100","99%"]
c5 <- SNHT.Critical.Values["100","95%"]
negC <- NegCurve(tSNHT)
abline(h=negC*c1,lty="dashed")
abline(h=negC*c5,lty="dotted")


dev.new()


# Figure 4
z<-ReadECAdata('RR_SOUID101991.txt')
d1900.2000<-z[which(z[,1]<'2001-01-01'),]
d1900.2000<-d1900.2000[which(d1900.2000[,1]>='1900-01-01'),]
y1900.2000w<-WetDayCount(d1900.2000,1,2,10)
y1900.2000 <- YearMeanFromDay(d1900.2000,1,2)
y1900.2000[,2]<-y1900.2000[,2]*10
original.parameters<-par()
par(mar=c(1,1,1,1),oma=c(2, 4, 1,4))
plot(y1900.2000[,1], y1900.2000[,2]*4, type = "s", axes = FALSE, ylim=c(0,1400), col="grey",xlab = "", ylab = "",main='HOMOGENEITY OF EUROPEAN CLIMATE SERIES')
mtext("Wet day count", side=2, line=3)
axis(side=4, col="grey",col.axis="grey")
par(new=T)
plot(y1900.2000w, type = "s", xlab = "", ylab = "",main='', axes = 1:2, ylim=c(80,320) )
mtext("Precipitation amount (mm)", side=4, line=3)
abline(v = c(1948), lty=2)
par(new=T)
x<-1:dim(y1900.2000)[1]
y<-y1900.2000[x,2]*4
y.loess <- loess(y ~ x, span=0.29, data.frame(x=x, y=y))
y.predict <- predict(y.loess, data.frame(x=x))
plot(y1900.2000[x,1],y.predict, type = "l",xlim=c(1906,2000),ylim=c(0,1400),axes = FALSE, col="grey",bty = "n",xlab = "", ylab = "")
par(new=T)
x<-1:dim(y1900.2000w)[1]
y<-y1900.2000w[x,2]
y.loess <- loess(y ~ x, span=0.29, data.frame(x=x, y=y))
y.predict <- predict(y.loess, data.frame(x=x))
plot(y1900.2000w[x,1],y.predict, type = "l",xlim=c(1906,2000),ylim=c(80,320),axes = FALSE, bty = "n",xlab = "", ylab = "")
par(original.parameters)

dev.new()


# Figure 5
y1900.2000w.1<-WetDayCount(d1900.2000,1,2,1)
y1900.2000w1<-WetDayCount(d1900.2000,1,2,10)
y1900.2000w10<-WetDayCount(d1900.2000,1,2,100)

wu<-ts(y1900.2000w.1[2:95,2],  start=c(y1900.2000w.1[1]), frequency=1)
BuishandRangeP.1<-BuishandRangeTest(wu)
wu<-ts(y1900.2000w1[2:95,2],  start=c(y1900.2000w1[1]), frequency=1)
BuishandRangeP1<-BuishandRangeTest(wu)
wu<-ts(y1900.2000w10[2:95,2],  start=c(y1900.2000w10[1]), frequency=1)
BuishandRangeP10<-BuishandRangeTest(wu)
par(mfrow=c(1,3))
plot(y1900.2000w1[2:95,1], BuishandRangeP.1, type = "l",ylim=c(-5,1),xlab = "", ylab = "Value statistic",main='Buishand range\nWet day count (threshold 0.1mm)')
plot(y1900.2000w1[2:95,1], BuishandRangeP1, type = "l",ylim=c(-5,1),xlab = "", ylab = "Value statistic",main='Buishand range\nWet day count (threshold 1mm)')
plot(y1900.2000w1[2:95,1], BuishandRangeP10, type = "l",ylim=c(-5,1),xlab = "", ylab = "Value statistic",main='Buishand range\nWet day count (threshold 10mm)')








