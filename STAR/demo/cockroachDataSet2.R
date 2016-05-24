require(STAR)
original <- par(ask = dev.interactive(orNone = TRUE))
oldpar <- par()

## load spontaneous data of 3 putative projection neurons
## simultaneously recorded from the cockroach (Periplaneta
## americana) antennal lobe
data(CAL2S)
## convert data into spikeTrain objects
CAL2S <- lapply(CAL2S,as.spikeTrain)
## look at the individual trains
## first the "raw" data
CAL2S[["neuron 1"]]
## next some summary information
summary(CAL2S[["neuron 1"]])
## next the renewal tests
renewalTestPlot(CAL2S[["neuron 1"]]);par(oldpar)
## It does not look too bad so let fit simple models
compModels(CAL2S[["neuron 1"]]);par(oldpar)
## the best one is the invgauss. Let's look at
## it in detail
isiHistFit(CAL2S[["neuron 1"]],"invgauss",xlim=c(0,0.5))
## refit the best model make a nicer looking
## QQ plot.
## WARNING the "MLE" fuctions take
## ISIs as arguments (not spikeTrains).
n1SDFig <- invgaussMLE(CAL2S[["neuron 1"]])
qqDuration(n1SDFig)
## The linear scale is perhaps not the best one
## so let's use the log scale
qqDuration(n1SDFig,log="xy")

## Do the same thing with neuron 2
CAL2S[["neuron 2"]]
## next some summary information
summary(CAL2S[["neuron 2"]])
## next the renewal tests
renewalTestPlot(CAL2S[["neuron 2"]]);par(oldpar)
compModels(CAL2S[["neuron 2"]]);par(oldpar)
## This time the invgauss is the best one
isiHistFit(CAL2S[["neuron 2"]],"invgauss",xlim=c(0,0.5))
n2SDFig <- invgaussMLE(CAL2S[["neuron 2"]])
qqDuration(n2SDFig)
qqDuration(n2SDFig,log="xy")

## Do the same thing with neuron 3
CAL2S[["neuron 3"]]
## next some summary information
summary(CAL2S[["neuron 3"]])
## next the renewal tests
renewalTestPlot(CAL2S[["neuron 3"]]);par(oldpar)
compModels(CAL2S[["neuron 3"]]);par(oldpar)
## This time the log logisitc is the best one
isiHistFit(CAL2S[["neuron 3"]],"llogis",xlim=c(0,0.5))
n3SDFll <- llogisMLE(CAL2S[["neuron 3"]])
qqDuration(n3SDFll)
qqDuration(n3SDFll,log="xy")

## build object for auto-correlation estimates
n1.1lt <- lockedTrain(CAL2S[[1]],,laglim=c(0,0.5))
n2.2lt <- lockedTrain(CAL2S[[2]],,laglim=c(0,0.5))
n3.3lt <- lockedTrain(CAL2S[[3]],,laglim=c(0,0.5))
## build object for cross-correlation estimates
n1.2lt <- lockedTrain(CAL2S[[1]],CAL2S[[2]],laglim=c(-0.5,0.5))
n1.3lt <- lockedTrain(CAL2S[[1]],CAL2S[[3]],laglim=c(-0.5,0.5))
n2.3lt <- lockedTrain(CAL2S[[2]],CAL2S[[3]],laglim=c(-0.5,0.5))
## display the whole thing on a matrix
layout(matrix(1:9,nrow=3))
## plot autoco on the diagonal, cross-raters below and
## cch above
## first column
hist(n1.1lt,bw=0.05,main="1.1",ylab="",type="l",col=2)
plot(n1.2lt,xlim=c(-0.1,0.1),main="1.2",ylab="",pch=".")
plot(n1.3lt,xlim=c(-0.1,0.1),main="1.3",ylab="",pch=".")
## second column
hist(n1.2lt,bw=0.05,main="1.2",ylab="",type="l",col=2)
hist(n2.2lt,bw=0.05,main="2.2",ylab="",type="l",col=2)
plot(n2.3lt,xlim=c(-0.1,0.1),main="2.3",ylab="",pch=".")
## third column
hist(n1.3lt,bw=0.05,main="1.3",ylab="",type="l",col=2)
hist(n2.3lt,bw=0.05,main="2.3",ylab="",type="l",col=2)
hist(n3.3lt,bw=0.05,main="3.3",ylab="",type="l",col=2)

## Look at the Citronellal responses
## Get the data
data(CAL2C)
## convert them into repeatedTrain objects
## The stimulus command is on between 5.87 s and 6.37s
CAL2C <- lapply(CAL2C,as.repeatedTrain)
## look at the individual raster plots
par(oldpar)
plot(CAL2C[["neuron 1"]],stimTimeCourse=c(5.87,6.37),main="N1")
plot(CAL2C[["neuron 2"]],stimTimeCourse=c(5.87,6.37),main="N2")
plot(CAL2C[["neuron 3"]],stimTimeCourse=c(5.87,6.37),main="N3")
## make psths with error bars on a single plot
layout(matrix(1:3,nrow=3))
psth(CAL2C[["neuron 1"]],stimTimeCourse=c(5.87,6.37),main="N1",
     breaks=25,colCI=2)
psth(CAL2C[["neuron 2"]],stimTimeCourse=c(5.87,6.37),main="N2",
     breaks=25,colCI=2)
psth(CAL2C[["neuron 3"]],stimTimeCourse=c(5.87,6.37),main="N3",
     breaks=25,colCI=2)
## make a finer psth for the neurons
psth(CAL2C[["neuron 1"]],stimTimeCourse=c(5.87,6.37),main="N1",
     breaks=c(0.25,0.05),colCI=2)
psth(CAL2C[["neuron 2"]],stimTimeCourse=c(5.87,6.37),main="N2",
     breaks=c(0.25,0.05),colCI=2)
psth(CAL2C[["neuron 3"]],stimTimeCourse=c(5.87,6.37),main="N3",
     breaks=c(0.25,0.05),colCI=2)

## get a summary for neuron 2 defining a response window between
## s6.5 and s7.5
summary(CAL2C[["neuron 2"]],c(6.5,7.5))
## check the p values for stationarity of the 3 neurons
sapply(CAL2C, function(l) summary(l)$globalPval)
## only neuron 4 seems to be stationary...
## make a cross-raster between neuron 1 and 3
plot(lockedTrain(CAL2C[["neuron 1"]],CAL2C[["neuron 2"]],laglim=0.01*c(-1,1)),stimTimeCourse=c(5.87,6.37),pch="*")
plot(lockedTrain(CAL2C[["neuron 1"]],CAL2C[["neuron 3"]],laglim=0.01*c(-1,1)),stimTimeCourse=c(5.87,6.37),pch="*")
plot(lockedTrain(CAL2C[["neuron 2"]],CAL2C[["neuron 3"]],laglim=0.01*c(-1,1)),stimTimeCourse=c(5.87,6.37),pch="*")

par(original)
