### R code from vignette source 'vignette.Rnw'

###################################################
### code chunk number 1: rarAnalysis1
###################################################

require("rtfbs")
# extract sequence and position weight matrix files from RTFBS package
exampleArchive <- system.file("extdata", "NRSF.zip", package="rtfbs")
unzip(exampleArchive, c("input.fas", "pwm.meme"))
# read the position weight matrices from meme text formatted file
pwm <- read.pwm("pwm.meme")
# delete the file since we now have the PWMs in memory
unlink("pwm.meme")


###################################################
### code chunk number 2: rarAnalysis2
###################################################
# read in the sequences
ms <- read.ms("input.fas")
# delete the file since we now have the sequences in memory
unlink("input.fas")


###################################################
### code chunk number 3: rarAnalysis3
###################################################

msSplit <- split.ms(ms, 500)


###################################################
### code chunk number 4: rarAnalysis4
###################################################

msGroups <- groupByGC.ms(msSplit, 4);


###################################################
### code chunk number 5: rarAnalyis5
###################################################

markovModels <- list()
for (i in 1:length(msGroups)) 
  markovModels[[i]] <- build.mm(msGroups[[i]], 3)


###################################################
### code chunk number 6: rarAnalysis6
###################################################

listOfScoresReal <- list()
for (i in 1:length(msGroups))
  listOfScoresReal[[i]] <- score.ms(msGroups[[i]], pwm, markovModels[[i]])


###################################################
### code chunk number 7: rarAnalysis8
###################################################

simulatedSeqs <- list()
for (i in 1:length(msGroups)) 
  simulatedSeqs[[i]] <- simulate.ms(markovModels[[i]], 500000)


###################################################
### code chunk number 8: rarAnalysis9
###################################################

listOfScoresSimulated <- list()
for (i in 1:length(msGroups))
  listOfScoresSimulated[[i]] <- score.ms(simulatedSeqs[[i]], pwm, markovModels[[i]])


###################################################
### code chunk number 9: rarAnalysis10
###################################################

fdrMap <- list()
for (i in 1:length(msGroups))
  fdrMap[[i]] <- calc.fdr(msGroups[[i]], listOfScoresReal[[i]], 
                          simulatedSeqs[[i]], listOfScoresSimulated[[i]])


###################################################
### code chunk number 10: Plot
###################################################

makeFdrPlot(fdrMap, col=rainbow(4))


###################################################
### code chunk number 11: rarAnalysis11
###################################################

bindingSites <- data.frame()
for (i in 1:length(msGroups)) {
  bindingSites <- rbind(bindingSites,
                        output.sites(listOfScoresReal[[i]],
                                     fdrScoreMap=fdrMap[[i]],
                                     fdrThreshold=0.1))
}


