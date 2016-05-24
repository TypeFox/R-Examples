require(STAR)
originalpar <- par(ask = dev.interactive(orNone = TRUE))
oldpar <- par()

## load spontaneous data of 4 putative projection neurons
## simultaneously recorded from the cockroach (Periplaneta
## americana) antennal lobe
data(CAL1S)
## convert data into spikeTrain objects
CAL1S <- lapply(CAL1S,as.spikeTrain)
## look at the individual trains
## first the "raw" data
CAL1S[["neuron 1"]]
## next some summary information
summary(CAL1S[["neuron 1"]])
## before generating the first "complicated" graph
## next the renewal tests
renewalTestPlot(CAL1S[["neuron 1"]]);par(oldpar)
## It does not look too bad so let fit simple models
compModels(CAL1S[["neuron 1"]]);par(oldpar)
## the best one is the invgauss. Let's look at
## it in detail
isiHistFit(CAL1S[["neuron 1"]],"invgauss",xlim=c(0,0.5));par(oldpar)
## refit the best model make a nicer looking
## QQ plot.
n1SDFig <- invgaussMLE(CAL1S[["neuron 1"]])
qqDuration(n1SDFig);par(oldpar)
## The linear scale is perhaps not the best one
## so let's use the log scale
qqDuration(n1SDFig,log="xy");par(oldpar)

## Do the same thing with neuron 2
CAL1S[["neuron 2"]]
## next some summary information
summary(CAL1S[["neuron 2"]])
## next the renewal tests
renewalTestPlot(CAL1S[["neuron 2"]]);par(oldpar)
compModels(CAL1S[["neuron 2"]]);par(oldpar)
## This time the log normal is the best one
isiHistFit(CAL1S[["neuron 2"]],"lnorm",xlim=c(0,0.5));par(oldpar)
n2SDFln <- lnormMLE(CAL1S[["neuron 2"]])
qqDuration(n2SDFln);par(oldpar)
qqDuration(n2SDFln,log="xy");par(oldpar)

## Do the same thing with neuron 3
CAL1S[["neuron 3"]]
## next some summary information
summary(CAL1S[["neuron 3"]])
## next the renewal tests
renewalTestPlot(CAL1S[["neuron 3"]]);par(oldpar)
compModels(CAL1S[["neuron 3"]]);par(oldpar)
## This time the log normal is the best one
isiHistFit(CAL1S[["neuron 3"]],"lnorm",xlim=c(0,0.5));par(oldpar)
n3SDFln <- lnormMLE(CAL1S[["neuron 3"]])
qqDuration(n3SDFln);par(oldpar)
qqDuration(n3SDFln,log="xy");par(oldpar)

## Do the same thing with neuron 4
CAL1S[["neuron 4"]]
## next some summary information
summary(CAL1S[["neuron 4"]])
## This neuron does not have enough spikes
## (50) for the renewal test...
compModels(CAL1S[["neuron 4"]]);par(oldpar)
## This time the gamma is the best one
isiHistFit(CAL1S[["neuron 4"]],"gamma",xlim=c(0,1.5));par(oldpar)
n4SDFga <- gammaMLE(CAL1S[["neuron 3"]])
qqDuration(n4SDFga);par(oldpar)
qqDuration(n4SDFga,log="xy");par(oldpar)

## build object for auto-correlation estimates
n1.1lt <- lockedTrain(CAL1S[[1]],,laglim=c(0,0.5))
n2.2lt <- lockedTrain(CAL1S[[2]],,laglim=c(0,0.5))
n3.3lt <- lockedTrain(CAL1S[[3]],,laglim=c(0,0.5))
n4.4lt <- lockedTrain(CAL1S[[4]],,laglim=c(0,0.5))
## build object for cross-correlation estimates
n1.2lt <- lockedTrain(CAL1S[[1]],CAL1S[[2]],laglim=c(-0.5,0.5))
n1.3lt <- lockedTrain(CAL1S[[1]],CAL1S[[3]],laglim=c(-0.5,0.5))
n1.4lt <- lockedTrain(CAL1S[[1]],CAL1S[[4]],laglim=c(-0.5,0.5))
n2.3lt <- lockedTrain(CAL1S[[2]],CAL1S[[3]],laglim=c(-0.5,0.5))
n2.4lt <- lockedTrain(CAL1S[[2]],CAL1S[[4]],laglim=c(-0.5,0.5))
n3.4lt <- lockedTrain(CAL1S[[3]],CAL1S[[4]],laglim=c(-0.5,0.5))
## display the whole thing on a matrix
layout(matrix(1:16,nrow=4))
## plot autoco on the diagonal, cross-raters below and
## cch above
## first column
hist(n1.1lt,bw=0.05,main="1.1",ylab="",type="l",col=2)
plot(n1.2lt,xlim=c(-0.1,0.1),main="1.2",ylab="")
plot(n1.3lt,xlim=c(-0.1,0.1),main="1.3",ylab="",pch=".")
plot(n1.4lt,xlim=c(-0.1,0.1),main="1.4",ylab="")
## second column
hist(n1.2lt,bw=0.05,main="1.2",ylab="",type="l",col=2)
hist(n2.2lt,bw=0.05,main="2.2",ylab="",type="l",col=2)
plot(n2.3lt,xlim=c(-0.1,0.1),main="2.3",ylab="")
plot(n2.4lt,xlim=c(-0.1,0.1),main="2.4",ylab="")
## third column
hist(n1.3lt,bw=0.05,main="1.3",ylab="",type="l",col=2)
hist(n2.3lt,bw=0.05,main="2.3",ylab="",type="l",col=2)
hist(n3.3lt,bw=0.05,main="3.3",ylab="",type="l",col=2)
plot(n3.4lt,xlim=c(-0.1,0.1),main="3.4",ylab="")
## fourth column
hist(n1.4lt,bw=0.05,main="1.4",ylab="",type="l",col=2)
hist(n2.4lt,bw=0.05,main="2.4",ylab="",type="l",col=2)
hist(n3.4lt,bw=0.05,main="3.4",ylab="",type="l",col=2)
hist(n4.4lt,bw=0.05,main="4.4",ylab="",type="l",col=2)

## Look at the Vanillin responses
## Get the data
data(CAL1V)
## convert them into repeatedTrain objects
## The stimulus command is on between 4.49 s and 4.99s
CAL1V <- lapply(CAL1V,as.repeatedTrain)
## look at the individual raster plots
par(oldpar)
plot(CAL1V[["neuron 1"]],stimTimeCourse=c(4.49,4.99),main="N1")
plot(CAL1V[["neuron 2"]],stimTimeCourse=c(4.49,4.99),main="N2")
plot(CAL1V[["neuron 3"]],stimTimeCourse=c(4.49,4.99),main="N3")
plot(CAL1V[["neuron 4"]],stimTimeCourse=c(4.49,4.99),main="N4")
## make psths with error bars on a single plot
layout(matrix(1:4,nrow=2))
psth(CAL1V[["neuron 1"]],stimTimeCourse=c(4.49,4.99),main="N1",
     breaks=20,colCI=2)
psth(CAL1V[["neuron 2"]],stimTimeCourse=c(4.49,4.99),main="N2",
     breaks=10,colCI=2)
psth(CAL1V[["neuron 3"]],stimTimeCourse=c(4.49,4.99),main="N3",
     breaks=10,colCI=2)
psth(CAL1V[["neuron 4"]],stimTimeCourse=c(4.49,4.99),main="N4",
     breaks=10,colCI=2)
## make a finer psth for neuron 1
par(oldpar)
psth(CAL1V[["neuron 1"]],stimTimeCourse=c(4.49,4.99),main="N1",
     breaks=c(0.25,0.05),colCI=2)
## get a summary for neuron 1 defining a response window between
## s5 and s6
summary(CAL1V[["neuron 1"]],c(5,6))
## check the p values for stationarity of the 4 neurons
sapply(CAL1V, function(l) summary(l)$globalPval)
## only neuron 4 seems to be stationary...
## make a cross-raster between neuron 1 and 3
plot(lockedTrain(CAL1V[["neuron 1"]],CAL1V[["neuron 3"]],laglim=0.01*c(-1,1)),stimTimeCourse=c(4.49,4.99),pch="*")
plot(lockedTrain(CAL1V[["neuron 1"]],CAL1V[["neuron 2"]],laglim=0.01*c(-1,1)),stimTimeCourse=c(4.49,4.99),pch="*")
plot(lockedTrain(CAL1V[["neuron 1"]],CAL1V[["neuron 4"]],laglim=0.01*c(-1,1)),stimTimeCourse=c(4.49,4.99),pch="*")

par(originalpar)
