### R code from vignette source 'Ch3.rnw'

###################################################
### code chunk number 1: setup
###################################################
source("GenericSettings.R")


###################################################
### code chunk number 2: Ch3.rnw:16-17
###################################################
data(Bike, package="MindOnStats")


###################################################
### code chunk number 3: Ch3.rnw:21-23
###################################################
summary(Bike)
summary(Bike$Type)


###################################################
### code chunk number 4: Ch3.rnw:28-30
###################################################
tapply(Bike$Time, Bike$Type, length)
tapply(Bike$Time, Bike$Type, length)/length(Bike$Type)


###################################################
### code chunk number 5: 1ABikePieChart
###################################################
pie(summary(Bike$Type))
title("Pie chart")


###################################################
### code chunk number 6: 1ABikeBarChart
###################################################
barplot(summary(Bike$Type))
title("Bar chart", ylab="Count", xlab="Type")


###################################################
### code chunk number 7: 1BBikePieChart
###################################################
pie(summary(Bike$Type), labels=paste(names(summary(Bike$Type)), round(100*summary(Bike$Type)/length(Bike$Type), 1)))
title("Pie chart")


###################################################
### code chunk number 8: 1BBikeBarChart
###################################################
barplot(summary(Bike$Type)/length(Bike$Type))
title("Bar chart", ylab="Percentage", xlab="Type")


###################################################
### code chunk number 9: Ch3.rnw:65-68
###################################################
tapply(Bike$Time, list(Bike$Gender, Bike$Type), length)
t(tapply(Bike$Time, list(Bike$Gender, Bike$Type), length))
addmargins(t(tapply(Bike$Time, list(Bike$Gender, Bike$Type), length)))


###################################################
### code chunk number 10: 2BikeStackBar
###################################################
barplot(tapply(Bike$Time, list(Bike$Gender, Bike$Type), length), legend=TRUE)
title("Stack bar chart of type, gender", ylab="Count", xlab="Type")


###################################################
### code chunk number 11: 2BikeClusterBar
###################################################
barplot(tapply(Bike$Time, list(Bike$Gender, Bike$Type), length), legend=TRUE, beside=TRUE)
title("Cluster bar chart of type, gender", ylab="Count", xlab="Type")


###################################################
### code chunk number 12: 3Bike3ClusterBar
###################################################
barplot(tapply(Bike$Time, list(paste(Bike$Gender, Bike$Direction), Bike$Type), length), legend=TRUE, beside=TRUE)
title("Chart of type, direction, gender", ylab="Count", xlab="Type")


###################################################
### code chunk number 13: Ch3.rnw:101-106
###################################################
Counts = tapply(Bike$Time, list(Bike$Gender, Bike$Type), length)
Counts
Divisor = matrix(tapply(Bike$Type, Bike$Type, length), nrow=2, ncol=4, byrow=T)
Divisor
round(t(100*Counts/Divisor), 1)


###################################################
### code chunk number 14: Ch3.rnw:114-116
###################################################
data(PennState1, package="MindOnStats")
str(PennState1)


###################################################
### code chunk number 15: Ch3.rnw:124-126
###################################################
data(HoldingBreath, package="MindOnStats")
str(HoldingBreath)


###################################################
### code chunk number 16: 9HBDotPlotWindows (eval = FALSE)
###################################################
## windows(7, 2.5)
## dotchart(HoldingBreath$Time, pch=20, xlab="Time (seconds)", lcolor=0)


###################################################
### code chunk number 17: 9HBDotPlotNoWindows
###################################################
dotchart(HoldingBreath$Time, pch=20, xlab="Time (seconds)", lcolor=0)


###################################################
### code chunk number 18: 9HBStripPlotJitterWindows (eval = FALSE)
###################################################
## windows(7, 2.5)
## stripchart(HoldingBreath$Time, method="jitter", pch=20, xlab="Time (seconds)")


###################################################
### code chunk number 19: 9HBStripPlotJitterNoWindows
###################################################
stripchart(HoldingBreath$Time, method="jitter", pch=20, xlab="Time (seconds)")


###################################################
### code chunk number 20: 9HBStripPlotStackWindows (eval = FALSE)
###################################################
## windows(7, 2.5)
## stripchart(round(HoldingBreath$Time,0), method="stack", pch=20, xlab="Time (seconds)")


###################################################
### code chunk number 21: 9HBStripPlotStackNoWindows
###################################################
stripchart(round(HoldingBreath$Time,0), method="stack", pch=20, xlab="Time (seconds)")


###################################################
### code chunk number 22: 10HBHist
###################################################
hist(HoldingBreath$Time, xlab="Time (seconds)")


###################################################
### code chunk number 23: HBStem
###################################################
stem(HoldingBreath$Time)


###################################################
### code chunk number 24: 12HBBoxPlotWindows (eval = FALSE)
###################################################
## windows(7, 5)
## boxplot(HoldingBreath$Time, horizontal=TRUE, xlab="Time (seconds)", main="Holding breath time")


###################################################
### code chunk number 25: 12HBBoxPlotNoWindows
###################################################
boxplot(HoldingBreath$Time, horizontal=TRUE, xlab="Time (seconds)", main="Holding breath time")


###################################################
### code chunk number 26: 16BikeIntPlot
###################################################
interaction.plot(x.factor=Bike$Type, trace.factor=Bike$Gender, response=Bike$Speed, ylab="Mean of Speed", xlab="Type", trace.label="Gender", main="Interaction plot for Speed")


###################################################
### code chunk number 27: 17HBScatterPlot
###################################################
attach(HoldingBreath)
plot(x=Age, y=Time, xlab="Age (years)", ylab="Time (seconds)", main="Holding Breath, Time vs Age")
detach(HoldingBreath)


###################################################
### code chunk number 28: Ch3.rnw:221-223
###################################################
data(Textbooks, package="MindOnStats")
str(Textbooks)


###################################################
### code chunk number 29: 18TextbooksScatterPlot
###################################################
plot(Price~Thickness, data=Textbooks, xlab="Thickness (mm)", ylab="Price (A$)", main="Textbook Price vs Thickness")


###################################################
### code chunk number 30: 20HBScatterPlotByGender
###################################################
attach(HoldingBreath)
plot(x=Age, y=Time, xlab="Age (years)", ylab="Time (seconds)", main="Holding Breath, Time vs Age", type="n")
points(x=Age[Gender=="m"], y=Time[Gender=="m"],  pch=21, col=4)
points(x=Age[Gender=="f"], y=Time[Gender=="f"],  pch=19, col=2)
legend("topright", title="Gender", legend=c("Male", "Female"), pch=c(21,19), col=c(4,2))
detach(HoldingBreath)


###################################################
### code chunk number 31: 21TextbooksScatterPlotByCover
###################################################
attach(Textbooks)
plot(Price~Thickness, xlab="Thickness (mm)", ylab="Price (A$)", main="Textbooks, Price vs Thickness", type="n")
points(Price[Coverstyle=="H"]~Thickness[Coverstyle=="H"],   pch=22, col="blue")
points(Price[Coverstyle=="S"]~Thickness[Coverstyle=="S"], pch=19, col="green")
legend("bottomright", title="Cover", legend=c("Hard", "Soft"), pch=c(22,19), col=c("blue","green"))
detach(Textbooks)


###################################################
### code chunk number 32: XDeathsTSPlot
###################################################
plot(ldeaths, ylim=c(0,4000), xlab="Year", ylab="Number of deaths", main="Deaths from lung diseases in the United Kingdom") 
lines(fdeaths, col="red")
lines(mdeaths, col="blue")
legend("topright", legend=c("Total", "Men", "Women"), col=c("black", "blue", "red"), lty=c(1,1,1))


