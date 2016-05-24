######################################################################
#  Demo of the code used in the book chapter
#  Hoehle, M. and A. Mazick, A. (2010) Aberration detection in R
#  illustrated by Danish mortality monitoring,  Book chapter in
#  T. Kass-Hout and X. Zhang (Eds.) Biosurveillance: A Health Protection
#  Priority, CRC Press.
#
#  The data read by csv files in the chapter are found as data("momo")
#  in the package. Courtesy to Statens Serum Institut for making
#  the mortality data public.
#
#  Author: Michael Hoehle
#  Date:   13 Oct 2009
######################################################################


#Load surveillance package
library("surveillance")

#Load Danish mortality data (see book chapter for CSV reading")
data("momo")

#Create a plot of the data as in Figure. 1 of the book chapter.
#Note: The year is determined by the ISO week, not the date
plot(momo[year(momo)>=2000,],ylab="No. of deaths",par.list=list(mar=c(4,2.2,2,1),cex.axis=1.5), type=observed ~ time | unit, col=c(gray(0.3),NA,NA),xaxis.tickFreq=list("%G"=atChange),xaxis.labelFormat="%G",xlab="time (weeks)")


par(mfrow=c(1,2),mar=c(4,4,2,1))
plot(momo,ylab="No. of deaths",xlab="time (weeks)",legend.opts=NULL, type=observed ~ time,col=c(gray(0.3),NA,NA),xaxis.tickFreq=list("%G"=atChange,"%m"=atChange),xaxis.labelFreq=list("%G"=atChange),xaxis.labelFormat="%G")
plot(momo[,"[0,1)"],xlab="time (weeks)",ylab="No. of deaths",legend.opts=NULL,col=c(gray(0.3),NA,NA),xaxis.tickFreq=list("%G"=atChange,"%m"=atChange),xaxis.labelFreq=list("%G"=atChange),xaxis.labelFormat="%G")
par(mfrow=c(1,1))

#Monitoring starts in week 40, 2007
phase2 <- which(epoch(momo) >= "2007-10-01")
s.far <- farrington(momo[,"[0,1)"], control=list(range=phase2,alpha=0.01,b=5,w=4,powertrans="none"))


cntrlFar <- s.far@control
upper.ptnone <-s.far@upperbound
cntrlFar$powertrans <- "2/3"
upper.pt23 <- farrington(momo[,"[0,1)"],control=cntrlFar)@upperbound
cntrlFar$powertrans <- "1/2"
upper.pt12 <- farrington(momo[,"[0,1)"],control=cntrlFar)@upperbound


## plot(s.far,ylab="No. of deaths",xlab="time (weeks)",main="")


ymax <- max(s.far@upperbound, upper.pt12, upper.pt23)*1.2
#par(mar=c(4,4,1,1))
plot(s.far,legend.opts=NULL,ylab="No. of deaths",main="",xlab="time (weeks)",ylim=c(0,ymax),col=c("darkgray",NA,gray(0.3)),lty=c(1,1,1),lwd=c(1,1,2),dx.upperbound=0,alarm.symbol=list(pch=24,col=1, cex=1))
lines(c(1:nrow(s.far)-0.5,nrow(s.far)+0.5),c(upper.pt12,upper.pt12[nrow(s.far)]),type="s",col="darkgray",lwd=2,lty=2)
lines(c(1:nrow(s.far)-0.5,nrow(s.far)+0.5),c(upper.pt23,upper.pt23[nrow(s.far)]),type="s",col=gray(0.1),lwd=2,lty=3)
legend(x="topright",c("none","1/2","2/3"),col=c(gray(0.3),"darkgray",gray(0.1)),lwd=2,lty=1:3,horiz=TRUE)
#legend(x="topright",c("none","1/2","2/3",expression(hat(mu)[t[0]])),col=c(gray(0.3),"darkgray",gray(0.1),1),lwd=c(2,2,2,3),lty=c(1:3,1),horiz=TRUE)

#Median of predictive distribution
lines(c(1:nrow(s.far)-0.5,nrow(s.far)+0.5),c(s.far@control$pd[,2],s.far@control$pd[nrow(s.far),2]),type="s",col=1,lwd=3)
text(nrow(s.far)+2,tail(observed(s.far),n=1),expression(hat(mu)[t[0]]))


alarmDates <- epoch(s.far[alarms(s.far) == 1,])


par(mar=c(4,4,2,2))
surv2 <- s.far
surv2@observed <- 0*surv2@observed
surv2@upperbound <- 0*surv2@observed

plot(surv2,ylim=c(-0.05,1),ylab="Quantile",xlab="time (weeks)",legend.opts=NULL,main="",dx.upperbound=0,alarm.symbol=list(pch=24,col=1, cex=1))
lines(surv2@control$pd[,1], type="S")
lines( c(1,nrow(surv2)+0.), rep( 1-s.far@control$alpha/2, 2),lty=2,col=1)


s.far.all <- farrington(momo, control=list(range=phase2,alpha=0.01,b=5,w=4))


## s.far.all <- farrington(momo, control=list(range=phase2,alpha=0.01,b=5,w=4))


## plot(s.far.all,type = alarm ~ time,xlab="time (weeks)")


par(mar=c(4,4,1,1))
plot(s.far.all,type = alarm ~ time,xlab="time (weeks)",main="",alarm.symbol=list(pch=24,col=1, cex=1.5),lvl=rep(1,nrow(s.far.all)))

#######################################################################
#Negative binomial GLM modelling using the population size as covariate
#######################################################################
phase1 <- which(year(momo) == 2002 & epochInYear(momo) == 40):(phase2[1]-1)
momo.df <- as.data.frame(momo)
m <- MASS::glm.nb( `observed.[75,85)` ~ 1 + epoch + sin(2*pi*epochInPeriod) + cos(2*pi*epochInPeriod) + `population.[75,85)`, data=momo.df[phase1,])
mu0 <- predict(m, newdata=momo.df[phase2,],type="response")


ci <- confint(m)


kappa <- 1.2
s.nb <- glrnb(momo[,"[75,85)"], control=list(range=phase2,alpha=1/m$theta,mu0=mu0,c.ARL=4.75,theta=log(kappa),ret="cases"))


alarmDates <- epoch(s.nb[alarms(s.nb) == 1,])


plot(s.nb,dx.upperbound=0,legend.opts=NULL,ylab="No. of deaths",main="",ylim=c(0,max(observed(s.nb))*1.1),xlab="time (weeks)",col=c("darkgray",NA,1),lwd=c(1,1,2),lty=c(1,1,1),alarm.symbol=list(pch=24,col=1, cex=1))
lines(mu0,lwd=2,col=1,lty=2)
lines(exp(log(mu0) + log(kappa)),col=1,lty=3,lwd=3)
legend(x=20,y=100,c(expression(mu[0,t]),expression(mu[1,t]),"NNBA"),col=c(1,1,1),lty=c(2,3,1),horiz=TRUE,bg="white",lwd=c(2,3,2))


set.seed(123)


######################################################################
# P(N_c <= 51|\tau=\infty) computation
######################################################################

#Number of simulations to perform. In book chapter this number is
#1000, but for the sake of a speedy illustration this is drastically
#reduced in this demonstration
nSims <- 10 #1000

######################################################################
# Simulate one run-length by first generating data from the negative
# binomial model and then applying the LR NegBin CUSUM to it
######################################################################

simone.TAleq65 <- function(sts, g) {
  observed(sts)[phase2,] <- rnbinom(length(mu0), mu=mu0, size=m$theta)
  one <- glrnb(sts, control=modifyList(control(s.nb), list(c.ARL=g)))
  return(any(alarms(one) > 0))
}

#Determine run-length using 1000 Monte Carlo samples
g.grid <- seq(1,8,by=0.5)
pMC <- sapply(g.grid, function(g) { 
  mean(replicate(nSims, simone.TAleq65(momo[,"[75,85)"],g))) 
})


#Density for comparison in the negative binomial distribution
dY <- function(y,mu,log=FALSE, alpha, ...) {
  dnbinom(y, mu=mu, size=1/alpha, log=log)
}
#nMax <- max(which( dY(0:1e4, mu=max(mu0),alpha=1/m$theta) >= 1e-20)) - 1


pMarkovChain <- sapply( g.grid, function(g) {
  TA <- LRCUSUM.runlength( mu=t(mu0), mu0=t(mu0), mu1=kappa*t(mu0), h=g, dfun = dY, n=rep(600,length(mu0)), alpha=1/m$theta)
  return(tail(TA$cdf,n=1))
})


par(mar=c(4,4,2,2))
matplot(g.grid, cbind(pMC,pMarkovChain),type="l",ylab=expression(P(T[A] <= 65 * "|" * tau * "=" * infinity)),xlab="g",col=1)
prob <- 0.1
lines(range(g.grid),rep(prob,2),lty=3,lwd=2)
axis(2,at=prob,las=1,cex.axis=0.7)
legend(x="topright",c("Monte Carlo","Markov chain"), lty=1:2,col=1)


m.01 <- MASS::glm.nb( `observed.[0,1)` ~ 1 + epoch + `population.[0,1)`+ sin(2*pi*epochInPeriod) + cos(2*pi*epochInPeriod), data=momo.df[phase1,])
mu0 <- predict(m.01, newdata=momo.df[phase2,],type="response")

#Correct for past outbreaks
#omega <- algo.farrington.assign.weights(residuals(m.01, type="deviance"))
#m.01.refit <- glm.nb( `observed.[0,1)` ~ 1 + epoch + `population.[0,1)`+ sin(2*pi*epochInPeriod) + cos(2*pi*epochInPeriod), data=momo.df[phase1,],weights=omega)
#mu0.refit <- predict(m.01.refit, newdata=momo.df[phase2,],type="response")

#Results from the previous Farrington method
mu0.far <- control(s.far)$pd[,2]


######################################################################
# Simulate one run-length by first generating data from the negative
# binomial model and then applying the LR NegBin CUSUM to it
######################################################################

simone.TAleq65.far <- function(sts, alpha, mu0, size) {
  observed(sts)[phase2,] <- rnbinom(length(mu0), mu=mu0, size=size)
  res <- farrington(sts, control=modifyList(control(s.far), list(alpha=alpha)))
  return(any(as.logical(alarms(res))))
}

#Determine run-length using 1000 Monte Carlo samples
res.far <- replicate(nSims, simone.TAleq65.far(momo[,"[0,1)"],alpha=0.01,mu0=mu0.far,size=m.01$theta))
(pTA65.far <- mean(res.far))

#Run CUSUM
kappa <- 1.2
s.nb.01 <- glrnb(momo[,"[0,1)"], control=list(range=phase2,alpha=1/m.01$theta,mu0=mu0.far,c.ARL=2.1,theta=log(kappa),ret="cases"))
alarmDates <- epoch(s.nb.01[alarms(s.nb.01) == 1,])

mu1 <- kappa*mu0.far
#Show as usual
plot(s.nb.01,dx.upperbound=0,legend.opts=NULL,ylab="No. of deaths",main="",xlab="time (weeks)",col=c("darkgray",NA,1),lwd=c(1,1,1),lty=c(1,1,1),ylim=c(0,max(s.nb.01@upperbound))*1.15,alarm.symbol=list(pch=24,col=1, cex=1))
lines(1:(nrow(s.far)+1)-0.5, c(mu0.far,tail(mu0.far,n=1)),lwd=3,col=1,lty=1,type="s")
lines(1:(nrow(s.far)+1)-0.5, c(mu1,tail(mu1,n=1)),col=1,lty=3,lwd=3,type="s")
legend(x="topright",c(expression(mu[0,t]),expression(mu[1,t]),"NNBA"),col=c(1,1,1),lty=c(1,3,1),horiz=TRUE,bg="white",lwd=c(3,3,1))

