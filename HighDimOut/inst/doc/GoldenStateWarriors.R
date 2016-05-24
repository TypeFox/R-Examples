## ----_1, echo=FALSE------------------------------------------------------
library(HighDimOut)
data(GoldenStatesWarriors)
summary(GoldenStatesWarriors)

## ----_2------------------------------------------------------------------
library(plyr)
data.scale <- t(aaply(.data = as.matrix(GoldenStatesWarriors[,-1]), .margins = 2, .fun = function(x) (x-mean(x, na.rm = T))/sd(x, na.rm = T)))
summary(data.scale)

data.scale[is.na(data.scale)] <- 0


## ----_3, warning=FALSE---------------------------------------------------
score.ABOD <- Func.ABOD(data = data.scale, basic = T)
score.trans.ABOD <- Func.trans(raw.score = score.ABOD, method = "ABOD")

## ----_4------------------------------------------------------------------
GoldenStatesWarriors$Name[order(score.trans.ABOD, decreasing = T)[1:5]]
GoldenStatesWarriors[order(score.trans.ABOD, decreasing = T)[1:5],]

## ----_5------------------------------------------------------------------
score.SOD <- Func.SOD(data = data.scale, k.nn = 10, k.sel = 5, alpha = .8)
score.trans.SOD <- Func.trans(raw.score = score.SOD, method = "SOD")

GoldenStatesWarriors$Name[order(score.trans.SOD, decreasing = T)[1:5]]
GoldenStatesWarriors[order(score.trans.SOD, decreasing = T)[1:5],]

