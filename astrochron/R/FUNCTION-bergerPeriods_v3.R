### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### bergerPeriods: determined the predicted orbital periods using Berger et al. (1992)
###         and (SRM: March 23, 2012; April 25, 2013; May 20, 2013; July 26, 2013)
###
###########################################################################

bergerPeriods <- function (age=0,genplot=T)
{

cat("\n----- COMPUTING PREDICTED PERIODS FOR OBLIQUITY AND PRECESSION -----\n")

# set up output array
etp <- rep(0,11*5)
dim(etp) <- c(11,5)

# assign ages
etp[1,1]=0
etp[2,1]=50
etp[3,1]=100
etp[4,1]=150
etp[5,1]=200
etp[6,1]=250
etp[7,1]=300
etp[8,1]=350
etp[9,1]=400
etp[10,1]=450
etp[11,1]=500

# assign o1
etp[1,2] = 54.0
etp[2,2] = 52.1
etp[3,2] = 50.2
etp[4,2] = 48.5
etp[5,2] = 46.7
etp[6,2] = 45.0
etp[7,2] = 42.9
etp[8,2] = 40.7
etp[9,2] = 38.7
etp[10,2] = 36.8
etp[11,2] = 35.0

# assign o2
etp[1,3] = 41.0
etp[2,3] = 39.9
etp[3,3] = 38.8
etp[4,3] = 37.7
etp[5,3] = 36.6
etp[6,3] = 35.6
etp[7,3] = 34.2
etp[8,3] = 32.9
etp[9,3] = 31.6
etp[10,3] = 30.3
etp[11,3] = 29.0

# assign p1
etp[1,4] = 23.0
etp[2,4] = 22.6
etp[3,4] = 22.3
etp[4,4] = 21.9
etp[5,4] = 21.5
etp[6,4] = 21.2
etp[7,4] = 20.7
etp[8,4] = 20.2
etp[9,4] = 19.7
etp[10,4] = 19.2
etp[11,4] = 18.7

# assign p2
etp[1,5] = 19.0 
etp[2,5] = 18.8
etp[3,5] = 18.5
etp[4,5] = 18.2
etp[5,5] = 18.0
etp[6,5] = 17.7
etp[7,5] = 17.4
etp[8,5] = 17.0
etp[9,5] = 16.7
etp[10,5] = 16.3
etp[11,5] = 16.0

etp <- data.frame(etp)

o1.lin <- approxfun(etp[,1],etp[,2],method="linear")
o2.lin <- approxfun(etp[,1],etp[,3],method="linear")
p1.lin <- approxfun(etp[,1],etp[,4],method="linear")
p2.lin <- approxfun(etp[,1],etp[,5],method="linear")

if(genplot)
 {
  par(mfrow=c(1,1))
  plot(1:500,p1.lin(1:500), type="l", col="red",ylim=c(15,54), xlab="Time (Ma)", ylab="Quasiperiod (ka)", main="Predicted Precession and Obliquity Periods of Berger et al. (1992)")
  points(etp[,1],etp[,4], cex=0.8,col="gray")
  points(age,p1.lin(age),col="red")
  text(age,p1.lin(age)+1,labels=c(p1.lin(age)),cex=0.8,col="red",font=2)
  lines(1:500,p2.lin(1:500), type="l", col="blue")
  points(age,p2.lin(age),col="blue")
  points(etp[,1],etp[,5], cex=0.8,col="gray")
  text(age,p2.lin(age)+1,labels=c(p2.lin(age)),cex=0.8,col="blue",font=2)
  lines(1:500,o1.lin(1:500), type="l", col="green")
  points(age,o1.lin(age),col="green")
  points(etp[,1],etp[,3], cex=0.8,col="gray")
  text(age,o1.lin(age)+1,labels=c(o1.lin(age)),cex=0.8,col="green",font=2)
  lines(1:500,o2.lin(1:500), type="l", col="black")
  points(age,o2.lin(age),col="black")
  points(etp[,1],etp[,2], cex=0.8,col="gray")
  text(age,o2.lin(age)+1,labels=c(o2.lin(age)),cex=0.8,col="black",font=2)
 }

cat(" * Predicted precession and obliquity periods (ka)=\n")
cat("     o1=", o1.lin(age),"\n")
cat("     o2=", o2.lin(age),"\n")
cat("     p1=", p1.lin(age),"\n")
cat("     p2=", p2.lin(age),"\n")

### END function bergerPeriods
}
