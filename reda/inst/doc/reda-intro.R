## ----setup---------------------------------------------------------------
library(reda) # attach package 
data(simuDat) # attach sample dataset

## ----data----------------------------------------------------------------
head(simuDat, 10)
str(simuDat)

## ----const---------------------------------------------------------------
constFit <- rateReg(Survr(ID, time, event) ~ group + x1, data = simuDat,
                    subset = ID %in% 1:50)
# brief summary
constFit # or explicitly call show(constFit)

## ----twoPieces-----------------------------------------------------------
# two pieces' constant rate function i.e. one internal knot
twoPiecesFit <- rateReg(Survr(ID, time, event) ~ group + x1, df = 2, 
                        data = simuDat, subset = ID %in% 1:50)
twoPiecesFit

## ----sixPieces-----------------------------------------------------------
piecesFit <- rateReg(Survr(ID, time, event) ~ group + x1, df = 2,
                     knots = seq(from = 28, to = 140, by = 28),
                     data = simuDat, subset = ID %in% 1:50)
piecesFit # note that df = 2 is neglected since knots are specified

## ----spline--------------------------------------------------------------
## df can be simply specified
splineFit <- rateReg(Survr(ID, time, event) ~ group + x1, df = 6,
                     degree = 3L, data = simuDat, subset = ID %in% 1:50)
## internal knots are set as 33% and 67% quantiles of time variable
splineFit 

## or internal knots are expicitly specified
splineFit <- rateReg(Survr(ID, time, event) ~ group + x1, df = 2,
                     degree = 3L, knots = c(56, 112),
                     data = simuDat, subset = ID %in% 1:50)
splineFit # note that df = 2 is neglected similarly

## ----summary-------------------------------------------------------------
summary(constFit)
summary(piecesFit, showCall = FALSE)
summary(splineFit, showCall = FALSE, showKnots = FALSE)

## ----est-----------------------------------------------------------------
## point estimates of covariate coefficients
coef(splineFit)
## confidence interval for covariate coefficients
confint(splineFit, level = 0.95) 
## estimated coefficients of baseline rate function
baseRate(splineFit)

## ------------------------------------------------------------------------
AIC(constFit, piecesFit, splineFit)
BIC(constFit, piecesFit, splineFit)

## ----sampleMcf-----------------------------------------------------------
## overall sample MCF
sampleMcf1 <- mcf(Survr(ID, time, event) ~ 1,
                  data = simuDat, subset = ID %in% 1:10)
## sample MCF for different groups
sampleMcf2 <- mcf(Survr(ID, time, event) ~ group,
                  data = simuDat, subset = ID %in% 1:10)

## ----plot:sampleMcf, fig.height = 5, fig.width = 7-----------------------
## plot overall sample MCF
plotMcf(sampleMcf1)
## plot MCF for different groups
plotMcf(sampleMcf2, mark.time = TRUE, 
        lty = c(1, 5), col = c("orange", "navy")) +
    ggplot2::xlab("Days") + ggplot2::theme_bw()

## ----piecesMcf, fig.height = 5, fig.width = 7----------------------------
piecesMcf <- mcf(piecesFit)
plotMcf(piecesMcf, conf.int = TRUE, col = "blueviolet") +
    ggplot2::xlab("Days") + ggplot2::theme_bw()

## ----splineMcf, fig.height = 5, fig.width = 7----------------------------
newDat <- data.frame(x1 = rep(0, 2), group = c("Treat", "Contr"))
estmcf <- mcf(splineFit, newdata = newDat, groupName = "Group", 
              groupLevels = c("Treatment", "Control"))
plotMcf(estmcf, conf.int = TRUE, col = c("royalblue", "red"), lty = c(1, 5)) +
    ggplot2::ggtitle("Control vs. Treatment") + ggplot2::xlab("Days") +
    ggplot2::theme_bw() 

