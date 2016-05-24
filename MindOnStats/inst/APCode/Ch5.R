### R code from vignette source 'Ch5.rnw'

###################################################
### code chunk number 1: setup
###################################################
source("GenericSettings.R")


###################################################
### code chunk number 2: Ch5.rnw:16-18
###################################################
data(PopSong, package="MindOnStats")
str(PopSong)


###################################################
### code chunk number 3: Ch5.rnw:22-24
###################################################
PopSong01 =PopSong[PopSong$Year==2001,]
str(PopSong01)


###################################################
### code chunk number 4: Ch5.rnw:27-31
###################################################
attach(PopSong01)
tapply(Genre, Genre, length)
length(Genre)
tapply(Genre, Genre, length)/length(Genre)


###################################################
### code chunk number 5: 1GenreBarPlot
###################################################
plot(Genre, xlab="Genre", ylab="Count", main="Counts of Songs by Genre for 2001")


###################################################
### code chunk number 6: Ch5.rnw:39-40
###################################################
detach(PopSong01)


###################################################
### code chunk number 7: Ch5.rnw:47-48
###################################################
aggregate(Position~Chart+Country+Author, data=PopSong, FUN=length)


###################################################
### code chunk number 8: Ch5.rnw:55-56
###################################################
pchisq(5.935, df=2, lower.tail=FALSE)


###################################################
### code chunk number 9: Ch5.rnw:60-61
###################################################
chisq.test(x=c(42, 49, 69), p=c(0.22, 0.4, 0.38))


###################################################
### code chunk number 10: Ch5.rnw:70-73
###################################################
attach(PopSong)
tapply(Position, list(Genre, Chart), length)
detach(PopSong)


###################################################
### code chunk number 11: Ch5.rnw:76-84
###################################################
attach(PopSong)
NewGenre = Genre
levels(NewGenre)
levels(NewGenre)[1]="AnHRock"
levels(NewGenre)[4]="AnHRock"
levels(NewGenre)
tapply(Position, list(NewGenre, Chart), length)
detach(PopSong)


###################################################
### code chunk number 12: Ch5.rnw:88-89
###################################################
chisq.test(NewGenre, PopSong$Chart)


###################################################
### code chunk number 13: Ch5.rnw:93-97
###################################################
BreastCancer = matrix(c(166,8340, 124, 7978), nrow=2, byrow=TRUE)
dimnames(BreastCancer)[[1]] = c("Hormones", "Placebo")
dimnames(BreastCancer)[[2]] = c("Yes", "No")
BreastCancer


###################################################
### code chunk number 14: Ch5.rnw:100-101
###################################################
chisq.test(BreastCancer)


