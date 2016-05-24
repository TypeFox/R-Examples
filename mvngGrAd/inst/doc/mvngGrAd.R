### R code from vignette source 'mvngGrAd.Rnw'

###################################################
### code chunk number 1: setWidth
###################################################
op <- options(); utils::str(op) 
options(width = 60)
options(prompt = "R> ")


###################################################
### code chunk number 2: loadPackage
###################################################
library(mvngGrAd)


###################################################
### code chunk number 3: sketchGrid1 (eval = FALSE)
###################################################
## 
## sketchGrid(i = 10,
##            j = 10,
##            rowLimit = 20,
##            colLimit = 20,
##            shapeCross =
##            list(1, ## down
##                 1, ## up
##                 1:4, ## left
##                 1:4),## right 
##            layers = NULL,
##            excludeCenter = TRUE)


###################################################
### code chunk number 4: sketchGrid2 (eval = FALSE)
###################################################
## 
## sketchGrid(i = 10,
##            j = 10,
##            rowLimit = 20,
##            colLimit = 20,
##            shapeCross =
##            list(2,
##                 2, ## up
##                 2:5, ## left
##                 2:5),## right 
##            layers = NULL,
##            excludeCenter = TRUE)


###################################################
### code chunk number 5: sketchGrid3 (eval = FALSE)
###################################################
##   sketchGrid(i = 10,
##            j = 10,
##            rowLimit = 20,
##            colLimit = 20,
##            shapeCross =
##            list(NULL, ## down
##                 NULL, ## up
##                 NULL, ## left
##                 NULL),## right 
##            layers = 1,
##            excludeCenter = TRUE)


###################################################
### code chunk number 6: sketchGrid4 (eval = FALSE)
###################################################
##   sketchGrid(i = 10,
##            j = 10,
##            rowLimit = 20,
##            colLimit = 20,
##            shapeCross =
##            list(NULL, ## down
##                 NULL, ## up
##                 NULL, ## left
##                 NULL),## right 
##            layers = 1:2,
##            excludeCenter = TRUE)


###################################################
### code chunk number 7: sketchGrid5 (eval = FALSE)
###################################################
## sketchGrid(i = 10,
##            j = 10,
##            rowLimit = 20,
##            colLimit = 20,
##            shapeCross =
##            list(1:4,
##                 1:4,
##                 1:4,
##                 1:4),
##            layers = c(1:4),
##            excludeCenter = TRUE)


###################################################
### code chunk number 8: sketchGrid5 (eval = FALSE)
###################################################
## sketchGrid(i = 10,
##            j = 10,
##            rowLimit = 20,
##            colLimit = 20,
##            shapeCross =
##            list(1:4,
##                 1:4,
##                 1:4,
##                 1:4),
##            layers = c(1:4),
##            excludeCenter = FALSE)


###################################################
### code chunk number 9: sketchGrid7 (eval = FALSE)
###################################################
## sketchGrid(i = 10,
##            j = 10,
##            rowLimit = 20,
##            colLimit = 20,
##            shapeCross =
##            list(2:4,
##                 2:4,
##                 2:4,
##                 2:4),
##            layers = c(2:4),
##            excludeCenter = TRUE)


###################################################
### code chunk number 10: sketchGrid8 (eval = FALSE)
###################################################
## 
## sketchGrid(i = 2,
##            j = 2,
##            rowLimit = 20,
##            colLimit = 20,
##            shapeCross =
##            list(2:4,
##                 2:4,
##                 2:4,
##                 2:4),
##            layers = c(2:4),
##            excludeCenter = TRUE)
## 


###################################################
### code chunk number 11: rowAndCol
###################################################
## row vector
rows <- rep(1 : 50, each = 50)

## column vector
cols <- rep(1 : 50, 50)


###################################################
### code chunk number 12: envError
###################################################
set.seed(13)

envError <-
  rep(c(seq(from = -12.5 , to = -0.5, by = 0.5),
        seq(from =   0.5 , to = 12.5, by = 0.5)),
      each = 50) + rnorm(2500)


###################################################
### code chunk number 13: scaleEnvError
###################################################
scaleFactE <-  (sd(envError) / sqrt(45))

envError <- 
  scale(envError, center = FALSE, 
        scale = scaleFactE)

envError <- as.vector(envError)


###################################################
### code chunk number 14: gEffects
###################################################
gEffects <- rnorm(2500,
                  mean = 0,
                  sd = 5)

scaleFactG <- (sd(gEffects)/sqrt(25))

gEffects <- 
  scale(gEffects,
        center = FALSE, 
        scale = scaleFactG)

gEffects <- as.vector(gEffects)


###################################################
### code chunk number 15: gValues
###################################################
## population mean
mu <- 30

gValues <- mu + gEffects


###################################################
### code chunk number 16: obsPisGplusE
###################################################
obsP <- gValues + envError


###################################################
### code chunk number 17: setNA
###################################################
obsP[3] <- NA


###################################################
### code chunk number 18: varPandVarG
###################################################
## phenotypic variance
varP <- var(obsP,na.rm = TRUE)

## genotypic variance
varG <- var(gEffects)


###################################################
### code chunk number 19: movingGrid1
###################################################
## creates object of class movG
resMG <- movingGrid(rows = rows,
                    columns = cols,
                    obsPhe = obsP,
                    shapeCross =
                    list(1:2,
                         1:2,
                         1:2,
                         1:2),
                    layers = 1)


###################################################
### code chunk number 20: summary
###################################################
summary(resMG)


###################################################
### code chunk number 21: entryData
###################################################
## only first 10
head(entryData(resMG))


###################################################
### code chunk number 22: fitted
###################################################
## only first 10
fitted(resMG)[1:10]


###################################################
### code chunk number 23: movingMean
###################################################
## only first 10
movingMean(resMG)[1:10]


###################################################
### code chunk number 24: residuals
###################################################
residuals(resMG)[1:10]


###################################################
### code chunk number 25: hsquareObs
###################################################
varG/varP


###################################################
### code chunk number 26: hsquareAdj
###################################################
## variance of the adj. phenot. values
varPadj <- var(fitted(resMG),na.rm = TRUE)

varG/varPadj


###################################################
### code chunk number 27: aToD
###################################################

 
layout(rbind(c(1,2),
             c(3,4)))
 
sketchGrid(i = 10,
           j = 10,
           rowLimit = 20,
           colLimit = 20,
           shapeCross =
           list(1, ## down
                1, ## up
                1:4, ## left
                1:4),## right 
           layers = NULL,
           excludeCenter = TRUE)

title(main = "[a]")

sketchGrid(i = 10,
           j = 10,
           rowLimit = 20,
           colLimit = 20,
           shapeCross =
           list(2, ## down
                2, ## up
                2:5, ## left
                2:5),## right 
           layers = NULL,
           excludeCenter = TRUE)

title(main = "[b]")

sketchGrid(i = 10,
           j = 10,
           rowLimit = 20,
           colLimit = 20,
           shapeCross =
           list(NULL, ## down
                NULL, ## up
                NULL, ## left
                NULL),## right 
           layers = 1,
           excludeCenter = TRUE)

title(main = "[c]")

sketchGrid(i = 10,
           j = 10,
           rowLimit = 20,
           colLimit = 20,
           shapeCross =
           list(NULL, ## down
                NULL, ## up
                NULL, ## left
                NULL),## right 
           layers = 1:2,
           excludeCenter = TRUE)

title(main = "[d]")



###################################################
### code chunk number 28: eToH
###################################################

layout(rbind(c(1,2),
             c(3,4)))

sketchGrid(i = 10,
           j = 10,
           rowLimit = 20,
           colLimit = 20,
           shapeCross =
           list(1:4,
                1:4,
                1:4,
                1:4),
           layers = c(1:4),
           excludeCenter = TRUE)

title(main = "[e]")


sketchGrid(i = 10,
           j = 10,
           rowLimit = 20,
           colLimit = 20,
           shapeCross =
           list(1:4,
                1:4,
                1:4,
                1:4),
           layers = c(1:4),
           excludeCenter = FALSE)

title(main = "[f]")

sketchGrid(i = 10,
           j = 10,
           rowLimit = 20,
           colLimit = 20,
           shapeCross =
           list(2:4,
                2:4,
                2:4,
                2:4),
           layers = c(2:4),
           excludeCenter = TRUE)

title(main = "[g]")

sketchGrid(i = 2,
           j = 2,
           rowLimit = 20,
           colLimit = 20,
           shapeCross =
           list(2:4,
                2:4,
                2:4,
                2:4),
           layers = c(2:4),
           excludeCenter = TRUE)

title(main = "[h]")


