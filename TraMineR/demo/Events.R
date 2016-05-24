require(grDevices); require(graphics)
oask <- devAskNewPage(dev.interactive(orNone = TRUE))

library(TraMineR)
data(actcal)
actcal.seq <- seqdef(actcal[,13:24])
transition <- seqetm(actcal.seq, method="transition")
transition[1,1:4] <- c("FullTime", "Decrease,PartTime",
		"Decrease,LowPartTime", "Stop")
transition[2,1:4] <- c("Increase,FullTime", "PartTime",
		"Decrease,LowPartTime", "Stop")
transition[3,1:4] <- c("Increase,FullTime", "Increase,PartTime",
		"LowPartTime", "Stop")
transition[4,1:4] <- c("Start,FullTime", "Start,PartTime",
		"Start,LowPartTime", "NoActivity")
print(transition)

actcal.seqe <- seqecreate(actcal.seq, tevent=transition)

fsubseq <- seqefsub(actcal.seqe, minSupport=100)
msubcount <- seqeapplysub(fsubseq, method="count")
#First lines...
msubcount[1:9, 6:7]
## Using time constraints
## Searching subsequences starting in summer (between June and September)
fsubseq <- seqefsub(actcal.seqe, minSupport=10,
		constraint=seqeconstraint(ageMin=6, ageMax=9))
fsubseq[1:10]
## Searching subsequences occurring in summer (between June and September)
fsubseq <- seqefsub(actcal.seqe, minSupport=10,
		constraint=seqeconstraint(ageMin=6, ageMax=9, ageMaxEnd=9))
fsubseq[1:10]
## Searching subsequences enclosed in a 6 months period
## and with a maximum gap of 2 months
fsubseq <- seqefsub(actcal.seqe, minSupport=10,
		constraint=seqeconstraint(maxGap=2, windowSize=6))
fsubseq[1:10]

## Looking for frequent subsequences
fsubseq <- seqefsub(actcal.seqe, pMinSupport=0.01)
## Frequences of 10 first subsequences
omar <- par(mar=c(5,4,4,2)+.1)
plot(fsubseq[1:10], col="cyan")
par(omar)
## looking for subsequence with FullTime
seqecontain(fsubseq, c("FullTime"))


## Looking for the discriminating subsequences for sex
sexsubseq <- seqecmpgroup(fsubseq, group=actcal$sex,
                    method="bonferroni")
## Plotting the ten most discriminating subsequences in 2 x 4 format
## frequencies
plot(sexsubseq[1:10],ptype="freq")
## residuals
plot(sexsubseq[1:10],ptype="resid")

devAskNewPage(oask)
