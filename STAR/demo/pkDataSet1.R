require(STAR)
originalpar <- par(ask = dev.interactive(orNone = TRUE))
oldpar <- par()

## load spontaneous data of 1 Purkinje cell
## recorded in cell attached mode from a cerebellar
## slice in control and bath applied bicuculline conditions
data(sPK)
## coerce data to spikeTrain objects
sPK <- lapply(sPK,as.spikeTrain)
## Get a summary of the ctl data
summary(sPK[["ctl"]])
## Look at the control train
## Don't show the rug plot for clarity
plot(sPK[["ctl"]],addRug=FALSE)
## Generate the renewal test plot taking into account
## the size of the data set (a lot of spikes!).
renewalTestPlot(sPK[["ctl"]],lag.max=250);par(oldpar)
## The data are stationary but NOT renewal. Look at the numerous
## nearly empty squares.

## Get a summary of the bicu data
summary(sPK[["bicu"]])
## Look at the control train
## Don't show the rug plot for clarity
plot(sPK[["bicu"]],addRug=FALSE)
## Generate the renewal test plot taking into account
## the size of the data set (a lot of spikes!).
renewalTestPlot(sPK[["bicu"]],lag.max=250);par(oldpar)
## This time the data are NOT stationary. This is seen clearly on a acf
## plot with very large lag.max
acf.spikeTrain(sPK[["bicu"]],lag.max=2000)


