### R code from vignette source 'pref.rnw'

###################################################
### code chunk number 1: pref.rnw:4-7
###################################################
library(lattice)
source("../../scripts/functions.R")
options(width=65)


###################################################
### code chunk number 2: PREF-base
###################################################
pref.tree.all <- read.csv("../../data/pref_trees.csv")


###################################################
### code chunk number 3: PREF-base
###################################################
names(pref.tree.all)
dim(pref.tree.all)


###################################################
### code chunk number 4: PREF-base
###################################################
show.cols.with.na(pref.tree.all)


###################################################
### code chunk number 5: PREF-base
###################################################
sapply(pref.tree.all, class)


###################################################
### code chunk number 6: PREF-base
###################################################
table(pref.tree.all$stratum)


###################################################
### code chunk number 7: PREF-base
###################################################
for (i in 1:4) 
  pref.tree.all[,i] <- factor(pref.tree.all[,i])


###################################################
### code chunk number 8: PREF-base
###################################################
mean(pref.tree.all$stratum)


###################################################
### code chunk number 9: PREF-base
###################################################
table(pref.tree.all$point)


###################################################
### code chunk number 10: PREF-base
###################################################
levels(pref.tree.all$point) <- c("1","3","5","7","9")


###################################################
### code chunk number 11: PREF-base
###################################################
table(pref.tree.all$point)


###################################################
### code chunk number 12: PREF-base
###################################################
pref.tree <- subset(pref.tree.all, dbh.cm > 25.4)
dim(pref.tree)


###################################################
### code chunk number 13: PREF-base
###################################################
levels(pref.tree$species)
pref.tree$species <- factor(pref.tree$species)
levels(pref.tree$species)


###################################################
### code chunk number 14: PREF-base
###################################################
pref.tree$vol.m3 <- pref.tree$vol.bf / 12 * 0.0283168466
pref.tree$tf.ha <- pref.tree$tf.ac * 2.47105381
pref.tree$vol.m3.ha <- pref.tree$vol.m3 * pref.tree$tf.ha
pref.tree$baf.m2.ha <- 
  pref.tree$baf.ft2.ac / 3.2808399^2 / 0.404685642


###################################################
### code chunk number 15: PREF-base
###################################################
pref.tree$weight <- 1/pref.tree$baf.ft2.ac


###################################################
### code chunk number 16: PREF-base
###################################################
pref.point <- with(pref.tree,
                   aggregate(x = list(ba.m2.ha = baf.m2.ha,
                                      vol.m3.ha = vol.m3.ha),
                             by = list(stratum = stratum,
                                       cluster = cluster,
                                       point = point),
                             FUN = sum))


###################################################
### code chunk number 17: PREF-base
###################################################
table(pref.point$stratum)


###################################################
### code chunk number 18: PREF-base
###################################################
design.point <- 
  expand.grid(cluster = levels(pref.tree.all$cluster),
              point = levels(pref.tree.all$point))


###################################################
### code chunk number 19: pref.rnw:290-291 (eval = FALSE)
###################################################
## str(design.point)


###################################################
### code chunk number 20: pref.rnw:293-294
###################################################
str(design.point, vec.len = 1)


###################################################
### code chunk number 21: PREF-base
###################################################
design.cluster <- 
  unique(pref.tree.all[, c("cluster","stratum")])


###################################################
### code chunk number 22: pref.rnw:311-312 (eval = FALSE)
###################################################
## str(design.cluster)


###################################################
### code chunk number 23: pref.rnw:314-315
###################################################
str(design.cluster, vec.len=1)


###################################################
### code chunk number 24: PREF-base
###################################################
dim(design.point)
dim(design.cluster)
test <- merge(x = design.point,
              y = design.cluster,
              all = TRUE)
dim(test)
head(test)


###################################################
### code chunk number 25: PREF-base
###################################################
design.point <- test
dim(design.point)
dim(pref.point)
test <- merge(x = design.point,
              y = pref.point,
              all = TRUE)
dim(test)
head(test)


###################################################
### code chunk number 26: PREF-base
###################################################
test$vol.m3.ha[is.na(test$vol.m3.ha)] <- 0
test$ba.m2.ha[is.na(test$ba.m2.ha)] <- 0
pref.point <- test
rm(test)


###################################################
### code chunk number 27: PREF-base
###################################################
pref.point$weight <- 1


###################################################
### code chunk number 28: PREF-pixel
###################################################
pref.pixel <- read.csv("../../data/pref_pixels.csv")
names(pref.pixel)
dim(pref.pixel)
pref.pixel[1:5, c(1:5,11)]
show.cols.with.na(pref.pixel)
sapply(pref.pixel, class)


###################################################
### code chunk number 29: PREF-pixel
###################################################
pref.subplot <- read.csv("../../data/pref_subplots.csv")
names(pref.subplot)
dim(pref.subplot)
pref.subplot[1:5, c(1:3, 5:6, 24)]
show.cols.with.na(pref.subplot)
sapply(pref.subplot, class)[1:4]
pref.subplot$plot <- factor(pref.subplot$plot)
pref.subplot$subplot <- factor(pref.subplot$subplot)


###################################################
### code chunk number 30: PREF-pixel
###################################################
dim(pref.point)
dim(pref.subplot[,c("plot","subplot","acndviC")])
pref.point.cov <- 
  merge(x = pref.point,
        y = pref.subplot[,c("plot","subplot","acndviC")],
        all.x = TRUE, all.y = FALSE,
        by.x = c("cluster","point"),
        by.y = c("plot","subplot"))
head(pref.point.cov)
dim(pref.point.cov)
show.cols.with.na(pref.point.cov)


###################################################
### code chunk number 31: ppwframe (eval = FALSE)
###################################################
## pref.pixel <- read.csv("../../data/pref_pixels.csv")
## print( wireframe( pref.pixel$elev ~ pref.pixel$utme + pref.pixel$utmn,
##           aspect = c(1, 0.5),
##           colorkey=FALSE,
##           pretty=TRUE,
##           shade=TRUE ) )


