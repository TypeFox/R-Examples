###################################################
### Chap13Start
###################################################
library(mistat)


###################################################
### DiceDesign
###################################################
library(lhs)

set.seed(123)

Des <- maximinLHS(n=14, k=7)

Des[, 1] <- Des[, 1] * (60-30) + 30 
Des[, 2] <- Des[, 2] * (0.02-0.005) + 0.005
Des[, 3] <- Des[, 3] * (0.01-0.002) + 0.002
Des[, 4] <- Des[, 4] * (5000-1000) + 1000
Des[, 5] <- Des[, 5] * (110000-90000) + 90000
Des[, 6] <- Des[, 6] * (296-290) + 290
Des[, 7] <- Des[, 7] * (360-340) + 340

Ps <- pistonSimulation(m=Des[,1], 
                       s=Des[,2], 
                       v0=Des[,3], 
                       k=Des[,4], 
                       p0=Des[,5], 
                       t=Des[,6], 
                       t0=Des[,7],
                       each=50, seed = 123)

Ps <- simulationGroup(Ps, 50)

aggregate(Ps[, !names(Ps) %in% "group"], by=Ps["group"], mean)


###################################################
### DiceLathyppiston2
###################################################
library(DiceEval)

data(LATHYPPISTON)

Dice <- modelFit(LATHYPPISTON[, !names(LATHYPPISTON) %in% "seconds"], 
                 LATHYPPISTON[,"seconds"], 
                 type = "Kriging", 
                 formula=~ ., 
                 control=list(trace=FALSE))

Dice$model

Dice <- modelFit(scale(x=LATHYPPISTON[, !names(LATHYPPISTON) %in% "seconds"]), 
                 LATHYPPISTON[,"seconds"], 
                 type = "Kriging", 
                 formula= ~ ., 
                 control=list(trace=FALSE))

Dice$model


###################################################
### PlotDiceLathyppiston2
###################################################
library(DiceView)

Dice <- km(design=LATHYPPISTON[, !names(LATHYPPISTON) %in% "seconds"], 
           response=LATHYPPISTON[,"seconds"], control=list(trace=FALSE))

sectionview(Dice, 
            center=colMeans(LATHYPPISTON[, !names(LATHYPPISTON) %in% "seconds"]), 
            conf_lev=c(0.5, 0.9, 0.95), 
            title="", col_sur="darkgrey", lwd=2, col_points="black", bg_blend=1,
            Xname=colnames(LATHYPPISTON[, !names(LATHYPPISTON) %in% "seconds"]))

layout(1)


###################################################
### Chap13End
###################################################
rm(LATHYPPISTON, Des, Ps, Dice)
detach(package:DiceView)
detach(package:DiceEval)
detach(package:DiceKriging)
detach(package:lhs)
