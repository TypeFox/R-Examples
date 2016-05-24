### R code from vignette source 'vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: vignette.Rnw:42-52
###################################################
options(width=65)
# load waterData package, assuming it has already been installed on the system
library(seawaveQ)
# load example data that comes with the package
data(swData)
# show first few rows of water-quality data for Missouri River at Omaha, Nebr.
head(qwMoRivOmaha)
# get a description of the data including definitions of the columns
# by viewing the help documentation
?qwMoRivOmaha


###################################################
### code chunk number 2: vignette.Rnw:61-63
###################################################
# scatter plot showing quantified, estimated, and censored  values
cenScatPlot(qwMoRivOmaha, pname="04035")


###################################################
### code chunk number 3: vignette.Rnw:70-85
###################################################
# scatter plot with many additional plotting arguments
# these options provide a plot closer to the plotting standards
# of the U.S. Geological Survey, however, these plots may not 
# meet all U.S. Geological Survey publication requirements
par(las=1, tcl=0.5)
cenScatPlot(qwMoRivOmaha, pname="04035", 
                       site="06610000 Missouri River at Omaha, Nebr.",
                       ylabel="Simazine concentration, in micrograms per liter",
                       legcex=0.7, qwcols=c("R", "P"), ylim=c(0,0.1), yaxs="i", 
                       cex.lab=0.9, cex.axis=0.9, xlim=c(as.Date("1996-01-01"), 
                       as.Date("2004-01-01")), xaxs="i", xaxt="n")
axdates <- c("1996-01-01", "1998-01-01", "2000-01-01", 
                       "2002-01-01", "2004-01-01")
axis(1, as.Date(axdates), 
                       labels=c("1996", "1998", "2000", "2002", "2004"), cex.axis=0.9)


###################################################
### code chunk number 4: vignette.Rnw:92-94
###################################################
# simple box plots of water-quality concentrations
rosBoxPlot(qwMoRivOmaha, qwcols=c("R", "P"))


###################################################
### code chunk number 5: vignette.Rnw:101-111
###################################################
# same boxplot function with many additional plotting arguments
rosBoxPlot(qwMoRivOmaha, site="06610000 Missouri River at Omaha, Nebr.",
                     log="y", yaxt="n", ylim=c(0.000001, 10), qwcols=c("R", "P"), 
                     ylab=c("Concentration, micrograms per liter"), col="skyblue1",
                     cex.axis=0.7, cex.sub=0.8, 
                     par(tcl=0.5, las=1, yaxs="i", mgp=c(3,0.5,0), mar=c(5,5,2,2), 
                     cex.main=0.9))
axis(2, at=c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10),
labels=c("0.000001", "0.00001", "0.0001", "0.001", "0.01",
"0.1", "1", "10"), cex.axis=0.7)


###################################################
### code chunk number 6: vignette.Rnw:118-124
###################################################
data(swData)
# show last few rows of water-quality data for Missouri River at Omaha, Nebr.
tail(cqwMoRivOmaha)
# get a description of the data including definitions of the columns
# by viewing the help documentation
?cqwMoRivOmaha


###################################################
### code chunk number 7: vignette.Rnw:133-138
###################################################
data(swData)
MoRivOmaha<-combineData(qwdat=qwMoRivOmaha, cqwdat=cqwMoRivOmaha,
qwcols=c("staid", "dates", "R", "P"))
# view combined data set
head(MoRivOmaha)


###################################################
### code chunk number 8: vignette.Rnw:146-165
###################################################
data(swData)

# associate continuous water-quality data with each sample
# combineData does this for you
modMoRivOmaha<-combineData(qwdat=qwMoRivOmaha, cqwdat=cqwMoRivOmaha)

# then fit model(s)
myfit1 <- fitswavecav(cdat=modMoRivOmaha, cavdat=cqwMoRivOmaha,
tanm="myfit1", pnames=c("04035", "04041"), yrstart=1995,
yrend=2003, tndbeg=1995, tndend=2003, iwcav=c("flowa30", "flowa1"),
dcol="dates", qwcols=c("R","P"))
myfit2 <- fitswavecav(cdat=modMoRivOmaha, cavdat=cqwMoRivOmaha,
tanm="myfit2", pnames=c("04035", "04041"), yrstart=1995,
yrend=2003, tndbeg=1995, tndend=2003, iwcav=c("seda30", "seda1"),
dcol="dates", qwcols=c("R","P"))
myfit3 <- fitswavecav(cdat=modMoRivOmaha, cavdat=cqwMoRivOmaha,
tanm="myfit3", pnames=c("04035", "04041"), yrstart=1995,
yrend=2003, tndbeg=1995, tndend=2003, iwcav=c("flowa30", "flowa1",
"seda30", "seda1"), dcol="dates", qwcols=c("R","P"))


###################################################
### code chunk number 9: vignette.Rnw:172-199
###################################################
# get the first element of the list for each model/constituent combination
# the data frame with information about each model/constituent combination
myfit1[[1]]
myfit2[[1]]
myfit3[[1]]

# get the second element of the list for each model/constituent combination
# the survival regression summary for each model/constituent combination
myfit1[[2]]
myfit2[[2]]
myfit3[[2]]

# get the first few lines of the third element of the list
head(myfit1[[3]])
head(myfit2[[3]])
head(myfit3[[3]])

# get the first few lines of the fourth element of the list
head(myfit1[[4]])
head(myfit2[[4]])
head(myfit3[[4]])

# get the summary of predicted concentrations
myfit1[[5]]
myfit2[[5]]
myfit3[[5]]



###################################################
### code chunk number 10: vignette.Rnw:206-211
###################################################

attributes(myfit1[[2]][[1]])
myfit1[[2]][[1]]$n
myfit1[[2]][[1]]$table



