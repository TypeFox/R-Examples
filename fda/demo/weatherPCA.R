#  Principal components analysis of the temperature and precipitatiion data

#  This file is intended to be used after the commands in files
#  weathersetup.R and weathersmooth.R have been executed.

#  Many other interesting ways of describing the data and plotting results
#  can be found in the file canadian-weather.R, set up in 2008 by 
#  Spencer Graves.

#  Last modified 17 November 2008

load("weatherfd")
load("weatherdata")

#  ---------  create fd objects for temp. and prec. ---------------

tempfdPar <- weatherfd$tempfdPar
precfdPar <- weatherfd$precfdPar

daytempfd <- tempfdPar$fd
dayprecfd <- precfdPar$fd

harmaccelLfd <- tempfdPar$Lfd

#  -----------------------------------------------------------------------
#               PCA of temperatures with varimax rotation
#  -----------------------------------------------------------------------

harmfdPar     <- fdPar(daytempfd, harmaccelLfd, 1e5)
daytemppcaobj <- pca.fd(daytempfd, nharm=4, harmfdPar)

daytemppcaobj <- varmx.pca.fd(daytemppcaobj)

#  plot harmonics

par(mfrow=c(1,1), pty="m")
plot.pca.fd(daytemppcaobj)

#  plot log eigenvalues

daytempeigvals <- daytemppcaobj[[2]]
par(ask=F)
plot(1:20, log10(daytempeigvals[1:20]), type="b",
     xlab="Eigenvalue Number", ylab="Log 10 Eigenvalue")
abline(lsfit(5:20, log10(daytempeigvals[5:20])), lty=2)

#  plot factor scores

harmscr <- daytemppcaobj[[3]]

plot(harmscr[,1], harmscr[,2],   xlab="Harmonic 3", ylab="Harmonic 4")
text(harmscr[,1], harmscr[,2]-4, weatherdata$station, cex=0.8, col=4)

#  -----------------------------------------------------------------------
#               PCA of precipitation with varimax rotation
#  -----------------------------------------------------------------------

harmfdPar     <- fdPar(dayprecfd, harmaccelLfd, 1e5)
dayprecpcaobj <- pca.fd(dayprecfd, nharm=3, harmfdPar)

dayprecpcaobj <- varmx.pca.fd(dayprecpcaobj)

#  plot harmonics

par(mfrow=c(1,1), pty="m")
plot.pca.fd(dayprecpcaobj)

#  plot log eigenvalues

daypreceigvals <- dayprecpcaobj[[2]]
par(ask=F)
plot(1:20, log10(daypreceigvals[1:20]), type="b",
     xlab="Eigenvalue Number", ylab="Log 10 Eigenvalue")
abline(lsfit(5:20, log10(daypreceigvals[5:20])), lty=2)

#  plot factor scores

harmscr <- dayprecpcaobj[[3]]

plot(harmscr[,1], harmscr[,2],   xlab="Harmonic 3", ylab="Harmonic 4")
text(harmscr[,1], harmscr[,2]-1, weatherdata$station, cex=0.8, col=4)



