### R code from vignette source 'trajVignette.Rnw'

###################################################
### code chunk number 1: trajVignette.Rnw:12-13
###################################################
options(width=80, digits=4, useFancyQuotes=FALSE, prompt=" ", continue=" ")


###################################################
### code chunk number 2: trajVignette.Rnw:31-33
###################################################
library(traj)
head(example.data$data)


###################################################
### code chunk number 3: trajVignette.Rnw:36-37
###################################################
head(example.data$time)


###################################################
### code chunk number 4: trajVignette.Rnw:73-75
###################################################
s1 = step1measures(example.data$data, example.data$time, ID = TRUE)
head(s1$measurments)


###################################################
### code chunk number 5: trajVignette.Rnw:80-82
###################################################
s2 = step2factors(s1)
head(s2$factors)


###################################################
### code chunk number 6: trajVignette.Rnw:87-88
###################################################
s3 = step3clusters(s2, nclusters = 4)


###################################################
### code chunk number 7: trajVignette.Rnw:93-95
###################################################
head(s3$clusters)
s3$clust.distr


###################################################
### code chunk number 8: trajVignette.Rnw:104-105
###################################################
plot(s3)


###################################################
### code chunk number 9: trajVignette.Rnw:115-116
###################################################
plotMeanTraj(s3)


###################################################
### code chunk number 10: trajVignette.Rnw:125-126
###################################################
plotMedTraj(s3)


###################################################
### code chunk number 11: trajVignette.Rnw:134-135
###################################################
plotBoxplotTraj(s3)


###################################################
### code chunk number 12: trajVignette.Rnw:143-144
###################################################
plotCombTraj(s3)


