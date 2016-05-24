### R code from vignette source 'mboost.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(prompt = "R> ", continue = "+  ", digits = 4)
set.seed(290875)


###################################################
### code chunk number 2: setup
###################################################
### make sure package mboost is available
if (!require("mboost"))
    install.packages("mboost", dependencies = TRUE)
library("mboost")
### package party is necessary for fitting trees
if (!require("party"))
    install.packages("party", dependencies = TRUE)
library("party")
### speed up things a little bit: if the model is already
### available load it into R
if (file.exists("model.Rda")) {
    load("model.Rda")
    mboost <- function(...) return(model)
    cvrisk <- function(...) return(cvm)
}


###################################################
### code chunk number 3: mboost-bodyfat-setup
###################################################
### set-up model formula
### names of features
data("bodyfat", package = "TH.data")
features <- names(bodyfat)[-2]
### set up model structure:
fml <- paste("bols(", features, ")", collapse = " + ") ### linear functions
fms <- paste("bbs(", features, ", center = TRUE, df = 1)",
             collapse = " + ")  ### smooth deviations from linearity
fmt <- "btree(hipcirc, waistcirc, tree_controls = ctree_control(maxdepth = 2, mincriterion = 0))" ### tree-based interaction
fm <- as.formula(paste("DEXfat", paste(fml, fms, fmt, sep = "+"), sep = "~"))


###################################################
### code chunk number 4: mboost-fm
###################################################
library("mboost")               ### attach package `mboost'
print(fm)                       ### model structure


###################################################
### code chunk number 5: mboost-bodyfat
###################################################
### fit model for conditional median of DEXfat
model <- mboost(fm,                           ### model structure
                data = bodyfat,               ### 71 observations
                family = QuantReg(tau = 0.5)) ### median regression


###################################################
### code chunk number 6: mboost-mstop
###################################################
model[1000]                     ### 900 more iterations


###################################################
### code chunk number 7: mboost-cvrisk
###################################################
### bootstrap for assessing the `optimal' number of boosting iterations
cvm <- cvrisk(model, grid = 1:100 * 10)
model[mstop(cvm)]               ### restrict model to optimal mstop(cvm) iterations


###################################################
### code chunk number 8: bodyfat-plot (eval = FALSE)
###################################################
## plot(cvm); plot(model)          ### depict out-of bag risk & selected components


###################################################
### code chunk number 9: mboost-bodyfat-plot-1
###################################################
### plot age and kneebreadth
cex <- 1.3
layout(matrix(c(1, 2, 1, 3), nr = 2))
par(mai = par("mai") * c(0.8, 1.1, 0.8, 0.8))
plot(cvm, cex.lab = cex)
mtext(text = "(A)", side = 3, 1, adj = 1)
plot(model, which = "bols(anthro3b", cex.lab = cex)
mtext(text = "(B)", side = 3, 1, adj = 1)
plot(model, which = "bbs(anthro3b",
     ylim = range(predict(model, which = "bols(anthro3b")),
     cex.lab = cex)
mtext(text = "(C)", side = 3, 1, adj = 1)


###################################################
### code chunk number 10: mboost-bodyfat-plot-2
###################################################
### plot interaction of hip and waist circumference
### first setup grid of hip and waist values
nd <- with(bodyfat,
           expand.grid(hipcirc = h <- seq(from = min(hipcirc),
                                          to = max(hipcirc),
                                          length = 100),
                  waistcirc = w <- seq(from = min(waistcirc),
                                       to = max(waistcirc),
                                       length = 100)))
### define colors for plot
col <-
c("#023FA5", "#1141A4", "#1A44A4", "#2146A4", "#2749A4", "#2C4BA4",
"#304DA4", "#3550A5", "#3852A5", "#3C54A6", "#4056A6", "#4359A7",
"#465BA7", "#495DA8", "#4C5FA9", "#4F61AA", "#5264AA", "#5566AB",
"#5868AC", "#5B6AAD", "#5D6CAE", "#606EAE", "#6270AF", "#6572B0",
"#6775B1", "#6A77B2", "#6C79B3", "#6F7BB4", "#717DB5", "#747FB6",
"#7681B6", "#7883B7", "#7B85B8", "#7D87B9", "#7F89BA", "#828BBB",
"#848DBC", "#868FBD", "#8891BE", "#8A93BE", "#8D94BF", "#8F96C0",
"#9198C1", "#939AC2", "#959CC3", "#979EC4", "#99A0C4", "#9BA1C5",
"#9DA3C6", "#9FA5C7", "#A1A7C8", "#A3A8C9", "#A5AAC9", "#A7ACCA",
"#A9AECB", "#ABAFCC", "#ACB1CC", "#AEB3CD", "#B0B4CE", "#B2B6CF",
"#B4B8CF", "#B5B9D0", "#B7BBD1", "#B9BCD2", "#BABED2", "#BCBFD3",
"#BEC1D4", "#BFC2D4", "#C1C4D5", "#C3C5D6", "#C4C7D6", "#C6C8D7",
"#C7C9D7", "#C9CBD8", "#CACCD9", "#CBCDD9", "#CDCFDA", "#CED0DA",
"#CFD1DB", "#D1D2DB", "#D2D3DC", "#D3D4DC", "#D4D6DD", "#D6D7DD",
"#D7D8DE", "#D8D9DE", "#D9DADF", "#DADBDF", "#DBDCDF", "#DCDCE0",
"#DDDDE0", "#DEDEE0", "#DEDFE1", "#DFDFE1", "#E0E0E1", "#E1E1E2",
"#E1E1E2", "#E2E2E2", "#E2E2E2", "#E2E2E2", "#E2E2E2", "#E2E2E2",
"#E2E2E2", "#E2E1E1", "#E2E0E1", "#E2E0E0", "#E1DFDF", "#E1DEDF",
"#E1DDDE", "#E1DCDD", "#E0DBDC", "#E0DADB", "#E0D9DA", "#DFD8D9",
"#DFD7D8", "#DFD6D7", "#DED5D6", "#DED3D5", "#DDD2D4", "#DDD1D3",
"#DDCFD2", "#DCCED0", "#DCCDCF", "#DBCBCE", "#DBCACD", "#DAC8CB",
"#DAC7CA", "#D9C5C8", "#D9C4C7", "#D8C2C6", "#D8C0C4", "#D7BFC3",
"#D7BDC1", "#D6BBC0", "#D5B9BE", "#D5B8BD", "#D4B6BB", "#D3B4B9",
"#D3B2B8", "#D2B0B6", "#D1AEB4", "#D1ADB3", "#D0ABB1", "#CFA9AF",
"#CEA7AE", "#CEA5AC", "#CDA3AA", "#CCA1A8", "#CB9FA7", "#CB9CA5",
"#CA9AA3", "#C998A1", "#C8969F", "#C7949D", "#C6929C", "#C5909A",
"#C48D98", "#C38B96", "#C38994", "#C28792", "#C18490", "#C0828E",
"#BF808C", "#BE7D8A", "#BD7B88", "#BB7986", "#BA7684", "#B97482",
"#B87180", "#B76F7E", "#B66C7C", "#B56A7A", "#B46777", "#B26575",
"#B16273", "#B05F71", "#AF5D6F", "#AE5A6D", "#AC576B", "#AB5569",
"#AA5266", "#A84F64", "#A74C62", "#A64960", "#A4475E", "#A3445B",
"#A24159", "#A03D57", "#9F3A55", "#9D3752", "#9C3450", "#9A304E",
"#992C4C", "#982949", "#962447", "#942045", "#931B42", "#911640",
"#900F3E", "#8E063B")
### use plot method to draw fitted values of the tree component only
print(plot(model, which = "btree", newdata = nd, col.regions = col,
     at = seq(from = -16, to = 16, length = 100)))
### save model for future use
save(model, cvm, file = "model.Rda")


###################################################
### code chunk number 11: mboost-predict
###################################################
### new data on a grid on range(anthro3b)
nd <- with(bodyfat, data.frame(anthro3b = seq(min(anthro3b), max(anthro3b),
                                              length = 100)))
### predictions for all base-learners of `anthro3b'
pr <- predict(model, which = "anthro3b", newdata = nd)
pr <- rowSums(pr)    ### aggregate linear and smooth effect


###################################################
### code chunk number 12: mboost-anthro3b
###################################################
plot(nd$anthro3b, pr, type = "l", xlab = "anthro3b",
     ylab = "f(anthro3b)")
lines(nd$anthro3b, predict(model, which = "bols(anthro3b", newdata = nd),
      type = "l", lty = "dashed")


