### R code from vignette source 'tutorial.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: tutorial.Rnw:19-22
###################################################
options(width=80)
options(continue=" ", prompt=" ")
RUNS <- 50 # for the simulations


###################################################
### code chunk number 2: installpackage (eval = FALSE)
###################################################
## install.packages(EloChoice) # install package (to be done once)


###################################################
### code chunk number 3: loadpackage
###################################################
library(EloChoice) # load package (every time you want to use the package)


###################################################
### code chunk number 4: tutorial.Rnw:56-59 (eval = FALSE)
###################################################
## xdata <- read.table("c:\\datafiles\\myfile.txt", sep="\t", header=T) # Windows
## xdata <- read.table("/Volumes/mydrive/myfile.txt", sep="\t", header=T) # Mac
## str(xdata)


###################################################
### code chunk number 5: tutorial.Rnw:64-74
###################################################
  library(xtable)
  winner <- c("ab", "cf", "ab", "dd", "ab")
  loser <- c("cf", "xs", "xs", "cf", "cf")
  rater <- c("A", "A", "A", "A", "B")
  date <- c("2010-01-01", "2010-01-01", "2010-01-01", "2010-01-01", "2010-01-04")
  time <- c("14:34:01", "14:34:08", "14:34:11", "14:34:15", "09:17:20")
  mytab <- data.frame(winner, loser, rater, date, time)
  colnames(mytab) <- c("preferred stimulus", "losing stimulus", "rater", "date", "time")
  cap <- "A possible data set layout. Note that R replaces space in column names with periods."
  print(xtable(mytab, align = , c("c", "c", "c", "c", "c", "c"), caption=cap,  label="tab:Ex"), include.rownames = FALSE, caption.placement="top")


###################################################
### code chunk number 6: createrandomdata
###################################################
set.seed(123)
xdata <- randompairs(nstim = 20, nint = 500, reverse = 0.1)
head(xdata)


###################################################
### code chunk number 7: tutorial.Rnw:101-104
###################################################
set.seed(123)
res <- elochoice(winner = xdata$winner, loser = xdata$loser, runs = 500)
summary(res)


###################################################
### code chunk number 8: tutorial.Rnw:113-114
###################################################
ratings(res, show="original", drawplot=FALSE)


###################################################
### code chunk number 9: tutorial.Rnw:118-119
###################################################
ratings(res, show="mean", drawplot=FALSE)


###################################################
### code chunk number 10: tutorial.Rnw:127-134 (eval = FALSE)
###################################################
## myratings <- ratings(res, show="mean", drawplot=FALSE)
## # Windows
## xdata <- write.table(myratings, "c:\\datafiles\\myratings.txt", sep="\t",
##                      header=T)
## # Mac
## xdata <- write.table(myratings, "/Volumes/mydrive/myratings.txt", sep="\t",
##                      header=T)


###################################################
### code chunk number 11: tutorial.Rnw:140-147 (eval = FALSE)
###################################################
## myratings <- ratings(res, show="all", drawplot=FALSE)
## # Windows
## xdata <- write.table(myratings, "c:\\datafiles\\myratings.txt", sep="\t",
##                      header=T, row.names=F)
## # Mac
## xdata <- write.table(myratings, "/Volumes/mydrive/myratings.txt", sep="\t",
##                      header=T, row.names=F)


###################################################
### code chunk number 12: fig1plot
###################################################
ratings(res, show=NULL, drawplot=TRUE)


###################################################
### code chunk number 13: fig1
###################################################
ratings(res, show=NULL, drawplot=TRUE)


###################################################
### code chunk number 14: tutorial.Rnw:181-190
###################################################
library(xtable)
pref <- c(1,1,0,0,1,0,1,0,0,0)
upset <- c("yes", "yes", "no", "no", "yes","no", "yes", "no","no","no" )
ratingdiff  <- c(200, 300, 100, 150, 200, 140, 280, 90, 150, 120)
ratingdiff2 <- c(90, 100, 300, 280, 120, 200, 140, 150, 200, 150)
mytab <- data.frame(pref, upset, ratingdiff, ratingdiff2)
colnames(mytab) <- c("higher rated = preferred", "upset", "rating difference 1", "rating difference 2")
cap <- "10 rating decisions that were either in accordance with the prediction or not. Two different rating differences are given to illustrate the weighted upset index. Note that the values are the same, just their assignment to different interactions is changed and consequently the column means are the same for both (173)."
print(xtable(mytab, digits = c(NA,0,NA,0, 0), align = , c("c", "c", "c", "c", "c"), caption=cap,  label="tab:upset"), include.rownames = FALSE, caption.placement="top")


###################################################
### code chunk number 15: tutorial.Rnw:197-199
###################################################
upsets <- reliability(res)
head(upsets)


###################################################
### code chunk number 16: tutorial.Rnw:206-208
###################################################
mean(upsets$upset)
mean(upsets$upset.wgt)


###################################################
### code chunk number 17: tutorial.Rnw:215-221
###################################################
set.seed(123)
xdata <- randompairs(nstim = 20, nint = 500, reverse = 0.3)
res <- elochoice(winner = xdata$winner, loser = xdata$loser, runs = 500)
upsets <- reliability(res)
mean(upsets$upset)
mean(upsets$upset.wgt)


###################################################
### code chunk number 18: tutorial.Rnw:233-238
###################################################
data(physical)
set.seed(123)
res <- elochoice(winner = physical$Winner, loser = physical$Loser, runs = 500)
summary(res)
ratings(res, show = "mean", drawplot = FALSE)


###################################################
### code chunk number 19: fig2plot
###################################################
ratings(res, show=NULL, drawplot=TRUE)


###################################################
### code chunk number 20: fig2
###################################################
ratings(res, show=NULL, drawplot=TRUE)


###################################################
### code chunk number 21: tutorial.Rnw:265-275
###################################################
# total of seven trials with two 'self-trials' (trials 6 and 7)
w <- c(letters[1:5], "a", "b"); l <- c(letters[2:6], "a", "b")
res <- elochoice(w, l)
ratings(res, drawplot=FALSE)
summary(res)
# total of five trials without 'self-trials'
w <- c(letters[1:5]); l <- c(letters[2:6])
res <- elochoice(w, l)
ratings(res, drawplot=FALSE)
summary(res)


