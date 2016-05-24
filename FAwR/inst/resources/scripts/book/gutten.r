### R code from vignette source 'gutten.rnw'

###################################################
### code chunk number 1: gutten.rnw:4-5
###################################################
source("../../scripts/functions.R")


###################################################
### code chunk number 2: gutten.rnw:18-22
###################################################
gutten <- read.csv("../../data/gutten.csv")
names(gutten) <- tolower(names(gutten))
names(gutten)[names(gutten)=="diameter"] <- "dbh.cm"
str(gutten)


###################################################
### code chunk number 3: gutten.rnw:25-26
###################################################
str(gutten, vec.len = 3)


###################################################
### code chunk number 4: gutten.rnw:31-34
###################################################
gutten$site <- factor(gutten$site)
gutten$location <- factor(gutten$location)
gutten$tree.ID <- with(gutten, interaction(location, tree))


###################################################
### code chunk number 5: gutten.rnw:43-47
###################################################
with(unique(gutten[,c("site","location","tree.ID")]),
     table(location))
with(unique(gutten[,c("site","location","tree.ID")]),
     table(location, site))


###################################################
### code chunk number 6: gutten.rnw:56-57
###################################################
show.cols.with.na(gutten)


###################################################
### code chunk number 7: gutten.rnw:68-69
###################################################
max(gutten$height[is.na(gutten$dbh.cm)])


###################################################
### code chunk number 8: gutten.rnw:74-75
###################################################
gutten <- gutten[!is.na(gutten$dbh.cm),]


