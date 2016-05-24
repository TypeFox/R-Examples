## ----lineWidthHook, eval=TRUE,echo = FALSE-------------------------------
library(knitr)
hook_output = knit_hooks$get('output')
knit_hooks$set(output = function(x, options) {
  # this hook is used only when the linewidth option is not NULL
  if (!is.null(n <- options$linewidth)) {
    x = knitr:::split_lines(x)
    # any lines wider than n should be wrapped
    if (any(nchar(x) > n)) x = strwrap(x, width = n)
    x = paste(x, collapse = '\n')
  }
  hook_output(x, options)
})

## ----installStMoMo, eval = FALSE-----------------------------------------
#  install.packages("StMoMo")

## ----loadStMoMo, eval=TRUE,echo = TRUE, message=FALSE, warning=FALSE-----
library(StMoMo)

## ----defineLClong, eval = FALSE------------------------------------------
#  constLC <- function(ax, bx, kt, b0x, gc, wxt, ages){
#    c1 <- mean(kt[1, ], na.rm = TRUE)
#    c2 <- sum(bx[, 1], na.rm = TRUE)
#    list(ax = ax + c1 * bx, bx = bx / c2, kt =  c2 * (kt - c1))
#  }
#  LC <- StMoMo(link = "logit", staticAgeFun = TRUE, periodAgeFun = "NP",
#               constFun = constLC)

## ----defineLC, eval=TRUE, cache=TRUE-------------------------------------
LC <- lc(link = "logit")

## ----defineCBDlong, eval = FALSE-----------------------------------------
#  f2 <- function(x, ages) x - mean(ages)
#  CBD <- StMoMo(link = "logit", staticAgeFun = FALSE,
#                periodAgeFun = c("1", f2))

## ----defineCBD, eval=TRUE, cache=TRUE------------------------------------
CBD <- cbd()

## ----defineAPC_RH_M7, eval=TRUE, cache=TRUE------------------------------
RH <- rh(link = "logit",  cohortAgeFun = "1")
APC <- apc(link = "logit")
M7 <- m7()

## ----definePlat, eval=TRUE, cache=TRUE-----------------------------------
f2 <- function(x, ages) mean(ages) - x
constPlat <- function(ax, bx, kt, b0x, gc, wxt, ages){
  nYears <- dim(wxt)[2]
  x <- ages
  t <- 1:nYears
  c <- (1 - tail(ages, 1)):(nYears - ages[1])
  xbar <- mean(x)
  #\sum g(c)=0, \sum cg(c)=0, \sum c^2g(c)=0
  phiReg <- lm(gc ~ 1 + c + I(c^2), na.action = na.omit)
  phi <- coef(phiReg)
  gc <- gc - phi[1] - phi[2] * c - phi[3] * c^2 
  kt[2, ] <- kt[2, ] + 2 * phi[3] * t
  kt[1, ] <- kt[1, ] + phi[2] * t + phi[3] * (t^2 - 2 * xbar * t)  
  ax <- ax + phi[1] - phi[2] * x + phi[3] * x^2
  #\sum kt[i, ] = 0
  ci <- rowMeans(kt, na.rm = TRUE)
  ax <- ax + ci[1] + ci[2] * (xbar - x) 
  kt[1, ] <- kt[1, ] - ci[1]
  kt[2, ] <- kt[2, ] - ci[2]
  list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
}
PLAT <- StMoMo(link = "logit", staticAgeFun = TRUE,
                  periodAgeFun = c("1", f2), cohortAgeFun = "1",
                  constFun = constPlat)

## ----showGNMformulaLC, eval=TRUE, cache=TRUE, dependson="defineLC", linewidth = 70----
LC$gnmFormula

## ----showGNMformulaCBD, eval=TRUE, cache=TRUE, dependson="defineCBD", linewidth = 70----
CBD$gnmFormula

## ----fit1, eval=TRUE, cache=TRUE, results='hide', dependson=c("defineLC", "defineAPC_RH_M7", "definePlat")----
Dxt <- EWMaleData$Dxt
Ext <- EWMaleData$Ext + 0.5 * EWMaleData$Dxt
ages <- EWMaleData$ages
years <- EWMaleData$years
ages.fit <- 55:89
wxt <- genWeightMat(ages = ages.fit, years = years, clip = 3)
LCfit <- fit(LC, Dxt = Dxt, Ext = Ext, ages = ages, years = years, 
             ages.fit = ages.fit, wxt = wxt)
APCfit <- fit(APC, Dxt = Dxt, Ext = Ext, ages = ages, years = years, 
             ages.fit = ages.fit, wxt = wxt)
CBDfit <- fit(CBD, Dxt = Dxt, Ext = Ext, ages = ages, years = years, 
             ages.fit = ages.fit, wxt = wxt)
M7fit <- fit(M7, Dxt = Dxt, Ext = Ext, ages = ages, years = years, 
             ages.fit = ages.fit, wxt = wxt)
PLATfit <- fit(PLAT, Dxt = Dxt, Ext = Ext, ages = ages, years = years, 
             ages.fit = ages.fit, wxt = wxt)

## ----fitRH, eval=TRUE, cache=TRUE, results='hide', dependson="defineAPC_RH_M7"----
RHfit <- fit(RH, Dxt = Dxt, Ext = Ext, ages = ages, years = years, 
             ages.fit = ages.fit, wxt = wxt, start.ax = LCfit$ax, 
             start.bx = LCfit$bx, start.kt = LCfit$kt)

## ----plotParam, eval=FALSE-----------------------------------------------
#  plot(LCfit, nCol = 3)
#  plot(CBDfit, parametricbx = FALSE)
#  plot(APCfit, parametricbx = FALSE, nCol = 3)

## ----plotLC, fig.width=8, fig.height = 2.8, out.width="16cm", echo=FALSE, cache=TRUE, dependson="fit1"----
plot(LCfit, nCol = 3)

## ----plotCBD, fig.width=8, fig.height = 4.2, out.width="11cm", echo=FALSE, cache=TRUE, dependson="fit1"----
plot(CBDfit, parametricbx = FALSE)

## ----plotAPC, fig.width=8, fig.height = 2.8, out.width="16cm", echo=FALSE, cache=TRUE, dependson="fit1"----
plot(APCfit, parametricbx = FALSE, nCol = 3)

## ----resAll, eval=TRUE, echo=FALSE, cache=TRUE, dependson=c("fit1","fitRH")----
LCres <- residuals(LCfit)
CBDres <- residuals(CBDfit)
APCres <- residuals(APCfit)
RHres <- residuals(RHfit)
M7res <- residuals(M7fit)
PLATres <- residuals(PLATfit)

## ----resLC_CBD, eval=FALSE-----------------------------------------------
#  LCres <- residuals(LCfit)
#  CBDres <- residuals(CBDfit)

## ----plotCBDresShow, eval=FALSE------------------------------------------
#  plot(CBDres, type = "colourmap", reslim = c(-3.5, 3.5))

## ----plotResShow, eval=FALSE---------------------------------------------
#  plot(LCres, type = "scatter", reslim = c(-3.5, 3.5))
#  plot(CBDres, type = "scatter", reslim = c(-3.5, 3.5))

## ----plotLCres, fig.width=4, fig.height=3,out.width="7.8cm", echo=FALSE, cache=TRUE, dependson="resAll"----
par(mar=c(4.5,4,1,1))
plot(LCres, type = "colourmap", reslim = c(-3.5,3.5))

## ----plotCBDres, fig.width=4, fig.height=3,out.width="7.8cm", echo=FALSE, cache=TRUE, dependson="resAll"----
par(mar=c(4.5,4,1,1))
plot(CBDres, type = "colourmap", reslim = c(-3.5,3.5))

## ----plotAPCres, fig.width=4, fig.height=3,out.width="7.8cm", echo=FALSE, cache=TRUE, dependson="resAll"----
par(mar=c(4.5,4,1,1))
plot(APCres, type = "colourmap", reslim = c(-3.5,3.5))

## ----plotRHres, fig.width=4, fig.height=3,out.width="7.8cm", echo=FALSE, cache=TRUE, dependson="resAll"----
par(mar=c(4.5,4,1,1))
plot(RHres, type = "colourmap", reslim = c(-3.5,3.5))

## ----plotM7res, fig.width=4, fig.height=3,out.width="7.8cm", echo=FALSE, cache=TRUE, dependson="resAll"----
par(mar=c(4.5,4,1,1))
plot(M7res, type = "colourmap", reslim = c(-3.5,3.5))

## ----plotPLATres, fig.width=4, fig.height=3,out.width="7.8cm", echo=FALSE, cache=TRUE, dependson="resAll"----
par(mar=c(4.5,4,1,1))
plot(PLATres, type = "colourmap", reslim = c(-3.5,3.5))

## ----plotLCresScatter, fig.width=10, fig.height=3,out.width="16cm", echo=FALSE, cache=TRUE, dependson="resAll"----
plot(LCres, type = "scatter", reslim = c(-3.5,3.5))

## ----plotCBDresScatter, fig.width=10, fig.height=3,out.width="16cm", echo=FALSE, cache=TRUE, dependson="resAll"----
plot(CBDres, type = "scatter", reslim = c(-3.5,3.5))

## ----BICCBD, cache=TRUE, dependson="fit1"--------------------------------
AIC(CBDfit)
BIC(CBDfit)

## ----BICtable, echo=FALSE, cache=TRUE, dependson="fit1", results='asis', message=FALSE, warning=FALSE----
library(xtable)
BICdf <-data.frame(list(LC = c(LCfit$npar[1], AIC(LCfit), BIC(LCfit)), 
                        CBD = c(CBDfit$npar[1], AIC(CBDfit), BIC(CBDfit)), 
                        APC = c(APCfit$npar[1], AIC(APCfit), BIC(APCfit)), 
                        RH = c(RHfit$npar[1], AIC(RHfit), BIC(RHfit)), 
                        M7 = c(M7fit$npar[1], AIC(M7fit), BIC(M7fit)), 
                        PLAT = c(PLATfit$npar[1], AIC(PLATfit), BIC(PLATfit))))
rownames(BICdf) <- c("Number of Parameters", "AIC", "BIC")
xtable(BICdf, dec = 1, caption = "Number of parameters, AIC and BIC values for different model fitted to the England and Wales males population for ages 55-89 and the period 1961-2011", center = "centering", file = "", floating = TRUE,label="tab:InfoCriteria",align="lcccccc", digits = 0)

## ----forAll, cache=TRUE, dependson=c("fit1","fitRH")---------------------
LCfor <- forecast(LCfit, h = 50)
CBDfor <- forecast(CBDfit, h = 50)
APCfor <- forecast(APCfit, h = 50, gc.order = c(1, 1, 0))
RHfor <- forecast(RHfit, h = 50, gc.order = c(1, 1, 0))
M7for <- forecast(M7fit, h = 50, gc.order = c(2, 0, 0))
PLATfor <- forecast(PLATfit, h = 50, gc.order = c(2, 0, 0))

## ----plotForecastPeriodIndexShow, eval=FALSE-----------------------------
#  plot(RHfor, only.kt = TRUE)
#  plot(M7for, only.kt = TRUE, nCol = 3)
#  plot(PLATfor, only.kt = TRUE)

## ----plotForecastCohortIndexShow, eval=FALSE-----------------------------
#  plot(APCfor, only.gc = TRUE)
#  plot(RHfor, only.gc = TRUE)
#  plot(M7for, only.gc = TRUE)
#  plot(PLATfor, only.gc = TRUE)

## ----simAll, cache=TRUE, dependson=c("fit1","fitRH")---------------------
set.seed(1234)
nsim <- 500
LCsim <- simulate(LCfit, nsim = nsim, h = 50)
CBDsim <- simulate(CBDfit, nsim = nsim, h = 50)
APCsim <- simulate(APCfit, nsim = nsim, h = 50, gc.order= c(1, 1, 0))
RHsim <- simulate(RHfit, nsim = nsim, h = 50, gc.order= c(1, 1, 0))
M7sim <- simulate(M7fit, nsim = nsim, h = 50, gc.order= c(2, 0, 0))
PLATsim <- simulate(PLATfit, nsim = nsim, h = 50, gc.order= c(2, 0, 0))

## ----plotForecastPeriodIndexRH, fig.width=5.2, fig.height=4,out.width="7.5cm", echo=FALSE, cache=TRUE, dependson="forAll"----
plot(RHfor, only.kt = TRUE, lwd = 2)

## ----plotForecastPeriodIndexPLAT, fig.width=10.4, fig.height=4,out.width="15cm", echo=FALSE, cache=TRUE, dependson="forAll"----
plot(PLATfor, only.kt = TRUE, lwd = 2)

## ----plotForecastPeriodIndexM7, fig.width=10.4, fig.height=2.8,out.width="22.5cm", echo=FALSE, cache=TRUE, dependson="forAll"----
plot(M7for, only.kt = TRUE, nCol = 3, lwd = 2)

## ----plotForecastCohortIndexAPC, fig.width=6, fig.height=4,out.width="8cm", echo=FALSE, cache=TRUE, dependson="forAll"----
plot(APCfor, only.gc = TRUE, lwd = 2)

## ----plotForecastCohortIndexRH, fig.width=6, fig.height=4,out.width="8cm", echo=FALSE, cache=TRUE, dependson="forAll"----
plot(RHfor, only.gc = TRUE, lwd = 2)

## ----plotForecastCohortIndexM7, fig.width=6, fig.height=4,out.width="8cm", echo=FALSE, cache=TRUE, dependson="forAll"----
plot(M7for, only.gc = TRUE, lwd = 2)

## ----plotForecastCohortIndexPLAT, fig.width=6, fig.height=4,out.width="8cm", echo=FALSE, cache=TRUE, dependson="forAll"----
plot(PLATfor, only.gc = TRUE, lwd = 2)

## ----plotSimulationRHShow, eval=FALSE------------------------------------
#  #Plot period index
#  par(mfrow = c(1, 3))
#  plot(RHfit$years, RHfit$kt[1, ],
#       xlim = range(RHfit$years, RHsim$kt.s$years),
#       ylim = range(RHfit$kt, RHsim$kt.s$sim[1, , 1:20]),
#       type = "l", xlab = "year", ylab = "kt", main = "Period index")
#  matlines(RHsim$kt.s$years, RHsim$kt.s$sim[1, , 1:20], type = "l", lty = 1)
#  #Plot cohort index
#  plot(RHfit$cohorts, RHfit$gc,
#       xlim = range(RHfit$cohorts, RHsim$gc.s$cohorts),
#       ylim = range(RHfit$gc, RHsim$gc.s$sim[, 1:20], na.rm = TRUE),
#       type = "l", xlab = "year", ylab = "kt",
#       main = "Cohort index (ARIMA(1,1,0) with drift)")
#  matlines(RHsim$gc.s$cohorts, RHsim$gc.s$sim[, 1:20], type = "l", lty = 1)
#  #Plot rates at age 65
#  qxt <- Dxt / Ext
#  plot(RHfit$years, qxt["65", ], xlim = range(RHfit$years, RHsim$years),
#       ylim = range(qxt["65", ], RHsim$rates["65", , 1:20]), type = "l",
#       xlab = "year", ylab = "rate", main = "Mortality rates at age 65")
#  matlines(RHsim$years, RHsim$rates["65", , 1:20], type = "l", lty = 1)

## ----plotSimulationRH, eval=TRUE,  fig.width=9, fig.height=3,out.width="16cm", results='hide', warning=FALSE, message=FALSE, cache=TRUE, echo=FALSE, dependson="simAll"----
par(mfrow=c(1, 3))
plot(RHfit$years, RHfit$kt[1, ],
     xlim = range(RHfit$years, RHsim$kt.s$years),
     ylim = range(RHfit$kt, RHsim$kt.s$sim[1, , 1:20]), 
     type = "l", xlab = "year", ylab = "kt", main = "Period index")
matlines(RHsim$kt.s$years, RHsim$kt.s$sim[1, , 1:20], type = "l", lty = 1)
#Plot cohort index
plot(RHfit$cohorts, RHfit$gc, 
     xlim = range(RHfit$cohorts, RHsim$gc.s$cohorts),
     ylim = range(RHfit$gc, RHsim$gc.s$sim[, 1:20], na.rm = TRUE),
     type = "l", xlab = "year", ylab = "kt",
     main = "Cohort index (ARIMA(1,1,0) with drift)")
matlines(RHsim$gc.s$cohorts, RHsim$gc.s$sim[, 1:20], type = "l", lty = 1)
#Plot rates at age 65
qxt <- Dxt / Ext
plot(RHfit$years, qxt["65", ], xlim = range(RHfit$years, RHsim$years), 
     ylim = range(qxt["65", ], RHsim$rates["65", , 1:20]), type = "l", 
     xlab = "year", ylab = "rate", main = "Mortality rates at age 65")
matlines(RHsim$years, RHsim$rates["65", , 1:20], type = "l", lty = 1)

## ----plotCBDFanShow, eval=FALSE------------------------------------------
#  library(fanplot)
#  probs = c(2.5, 10, 25, 50, 75, 90, 97.5)
#  qxt <- Dxt / Ext
#  #age 65
#  plot(CBDfit$years, qxt["65", ], xlim = c(1960, 2061), ylim = c(0.0025, 0.2),
#       xlab = "year", ylab = "mortality rate (log scale)",
#       pch = 20, log = "y")
#  fan(t(CBDsim$rates["65", , ]), start = 2012, probs = probs, n.fan = 4,
#      fan.col = colorRampPalette(c("black", "white")), ln = NULL)
#  #age 75
#  points(CBDfit$years, qxt["75", ], pch = 20)
#  fan(t(CBDsim$rates["75", , ]), start = 2012, probs = probs, n.fan = 4,
#      fan.col = colorRampPalette(c("red", "white")), ln = NULL)
#  #age 75
#  points(CBDfit$years, qxt["85", ], pch = 20)
#  fan(t(CBDsim$rates["85", , ]), start = 2012, probs = probs, n.fan = 4,
#      fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
#  #labels
#  text(1965, qxt[c("65", "75", "85"), "1990"],
#       labels = c("x = 65", "x = 75", "x = 85"))

## ----plotLCfan, fig.width=4, fig.height=5,out.width="4.5cm", echo=FALSE, warning=FALSE, cache=TRUE, dependson="simAll"----
library(fanplot)
probs = c(2.5, 10, 25, 50, 75, 90, 97.5)
qxt <- Dxt/Ext 
par(mar=c(4.5,4,1,1))
#age 65
plot(LCfit$years, qxt["65", ], xlim = c(1960, 2061), ylim = c(0.0025,0.2),
 xlab = "year", ylab = "mortality rate (log scale)", pch = 20, log = "y")
fan(t(LCsim$rates["65", , ]), start = 2012,
probs = probs, n.fan = 4,
fan.col = colorRampPalette(c("black", "white")), ln = NULL)
#age 75
points(LCfit$years, qxt["75", ], pch = 20)
fan(t(LCsim$rates["75", , ]), start = 2012,
probs = probs, n.fan = 4,
fan.col = colorRampPalette(c("red", "white")), ln = NULL)
#age 75
points(LCfit$years, qxt["85",], pch = 20)
fan(t(LCsim$rates["85", , ]), start = 2012,
probs = probs, n.fan = 4,
fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
#labels
text(1965, qxt[c("65","75", "85"),"1990"], labels=c("x=65", "x=75", "x=85"))

## ----plotCBDfan, fig.width=4, fig.height=5,out.width="4.5cm", echo=FALSE, cache=TRUE, dependson="plotLCfan"----
par(mar=c(4.5,4,1,1))
#age 65
plot(CBDfit$years, qxt["65", ], xlim = c(1960, 2061), ylim = c(0.0025,0.2),
     xlab = "year", ylab = "mortality rate (log scale)", pch = 20, log = "y")
fan(t(CBDsim$rates["65", , ]), start = 2012,
    probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("black", "white")), ln = NULL)
#age 75
points(CBDfit$years, qxt["75", ], pch = 20)
fan(t(CBDsim$rates["75", , ]), start = 2012,
    probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("red", "white")), ln = NULL)
#age 75
points(CBDfit$years, qxt["85",], pch = 20)
fan(t(CBDsim$rates["85", , ]), start = 2012,
    probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
#labels
text(1965, qxt[c("65","75", "85"),"1990"], labels=c("x=65", "x=75", "x=85"))

## ----plotAPCfan, fig.width=4, fig.height=5,out.width="4.5cm", echo=FALSE, cache=TRUE, dependson="plotLCfan"----
par(mar=c(4.5,4,1,1))
#age 65
plot(APCfit$years, qxt["65", ], xlim = c(1960, 2061), ylim = c(0.0025,0.2),
 xlab = "year", ylab = "mortality rate (log scale)", pch = 20, log = "y")
fan(t(APCsim$rates["65", , ]), start = 2012,
probs = probs, n.fan = 4,
fan.col = colorRampPalette(c("black", "white")), ln = NULL)
#age 75
points(APCfit$years, qxt["75", ], pch = 20)
fan(t(APCsim$rates["75", , ]), start = 2012,
probs = probs, n.fan = 4,
fan.col = colorRampPalette(c("red", "white")), ln = NULL)
#age 75
points(APCfit$years, qxt["85",], pch = 20)
fan(t(APCsim$rates["85", , ]), start = 2012,
probs = probs, n.fan = 4,
fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
#labels
text(1965, qxt[c("65","75", "85"),"1990"], labels=c("x=65", "x=75", "x=85"))

## ----plotRHfan, fig.width=4, fig.height=5,out.width="4.5cm", echo=FALSE, cache=TRUE, dependson="plotLCfan"----
par(mar=c(4.5,4,1,1))
#age 65
plot(RHfit$years, qxt["65", ], xlim = c(1960, 2061), ylim = c(0.0025,0.2),
 xlab = "year", ylab = "mortality rate (log scale)", pch = 20, log = "y")
fan(t(RHsim$rates["65", , ]), start = 2012,
probs = probs, n.fan = 4,
fan.col = colorRampPalette(c("black", "white")), ln = NULL)
#age 75
points(RHfit$years, qxt["75", ], pch = 20)
fan(t(RHsim$rates["75", , ]), start = 2012,
probs = probs, n.fan = 4,
fan.col = colorRampPalette(c("red", "white")), ln = NULL)
#age 85
points(RHfit$years, qxt["85",], pch = 20)
fan(t(RHsim$rates["85", , ]), start = 2012,
probs = probs, n.fan = 4,
fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
#labels
text(1965, qxt[c("65","75", "85"),"1990"], labels=c("x=65", "x=75", "x=85"))

## ----plotM7fan, fig.width=4, fig.height=5,out.width="4.5cm", echo=FALSE, cache=TRUE, dependson="plotLCfan"----
par(mar=c(4.5,4,1,1))
#age 65
plot(M7fit$years, qxt["65", ], xlim = c(1960, 2061), ylim = c(0.0025,0.2),
 xlab = "year", ylab = "mortality rate (log scale)", pch = 20, log = "y")
fan(t(M7sim$rates["65", , ]), start = 2012,
probs = probs, n.fan = 4,
fan.col = colorRampPalette(c("black", "white")), ln = NULL)
#age 75
points(M7fit$years, qxt["75", ], pch = 20)
fan(t(M7sim$rates["75", , ]), start = 2012,
probs = probs, n.fan = 4,
fan.col = colorRampPalette(c("red", "white")), ln = NULL)
#age 85
points(M7fit$years, qxt["85",], pch = 20)
fan(t(M7sim$rates["85", , ]), start = 2012,
probs = probs, n.fan = 4,
fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
#labels
text(1965, qxt[c("65","75", "85"),"1990"], labels=c("x=65", "x=75", "x=85"))

## ----plotPLATfan, fig.width=4, fig.height=5,out.width="4.5cm", echo=FALSE, cache=TRUE, dependson="plotLCfan"----
par(mar=c(4.5,4,1,1))
#age 65
plot(PLATfit$years, qxt["65", ], xlim = c(1960, 2061), ylim = c(0.0025,0.2),
 xlab = "year", ylab = "mortality rate (log scale)", pch = 20, log = "y")
fan(t(PLATsim$rates["65", , ]), start = 2012,
probs = probs, n.fan = 4,
fan.col = colorRampPalette(c("black", "white")), ln = NULL)
#age 75
points(PLATfit$years, qxt["75", ], pch = 20)
fan(t(PLATsim$rates["75", , ]), start = 2012,
probs = probs, n.fan = 4,
fan.col = colorRampPalette(c("red", "white")), ln = NULL)
#age 85
points(PLATfit$years, qxt["85",], pch = 20)
fan(t(PLATsim$rates["85", , ]), start = 2012,
probs = probs, n.fan = 4,
fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
#labels
text(1965, qxt[c("65","75", "85"),"1990"], labels=c("x=65", "x=75", "x=85"))

## ----loadNZ, eval=FALSE--------------------------------------------------
#    library(demography)
#    NZdata <- hmd.mx(country = "NZL_NP", username = username,
#                     password = password)

## ----fitNZ, eval=FALSE---------------------------------------------------
#    Ext_NZ <- NZdata$pop$male
#    Dxt_NZ <- NZdata$rate$male * Ext_NZ
#    LCfit_NZ <- fit(lc(), Dxt = Dxt_NZ, Ext = Ext_NZ, ages = NZdata$age,
#               years = NZdata$year, ages.fit = 0:89, years.fit = 1985:2008)

## ----bootLCNZ, eval=FALSE------------------------------------------------
#  LCboot_NZ <- bootstrap(LCfit_NZ, nBoot = 5000, type = "semiparametric")

## ----PlotbootLCNZshow, eval=FALSE----------------------------------------
#  plot(LCboot_NZ)

## ----simPULCNZ, eval=FALSE-----------------------------------------------
#  LCsimPU_NZ <- simulate(LCboot_NZ,  h = 24)

## ----simLCNZ, eval=FALSE-------------------------------------------------
#  LCfor_NZ <- forecast(LCfit_NZ, h = 24)
#  LCsim_NZ <- simulate(LCfit_NZ, nsim = 5000, h = 24)

## ----plotLCPredshow, eval=FALSE------------------------------------------
#  #Observed, fitted and central forecasts
#  mxt <- LCfit_NZ$Dxt / LCfit_NZ$Ext
#  mxtHat <- fitted(LCfit_NZ, type = "rates")
#  mxtCentral <- LCfor_NZ$rates
#  #95% Prediction intervals without parameter uncertainty
#  mxtPred2.5 <- apply(LCsim_NZ$rates, c(1, 2), quantile, probs = 0.025)
#  mxtPred97.5 <- apply(LCsim_NZ$rates, c(1, 2), quantile, probs = 0.975)
#  #95% intervals with parameter uncertainty (in sample, and predictions)
#  mxtHatPU2.5 <- apply(LCsimPU_NZ$fitted, c(1, 2), quantile, probs = 0.025)
#  mxtHatPU97.5 <- apply(LCsimPU_NZ$fitted, c(1, 2), quantile, probs = 0.975)
#  mxtPredPU2.5 <- apply(LCsimPU_NZ$rates, c(1, 2), quantile, probs = 0.025)
#  mxtPredPU97.5 <- apply(LCsimPU_NZ$rates, c(1, 2), quantile, probs = 0.975)
#  #Plot
#  x <- c("40", "60", "80")
#  matplot(LCfit_NZ$years, t(mxt[x, ]),
#          xlim = range(LCfit_NZ$years, LCfor_NZ$years),
#          ylim = range(mxtHatPU97.5[x, ], mxtPredPU2.5[x, ], mxt[x, ]),
#          type = "p", xlab = "years", ylab = "mortality rates (log scale)",
#          log = "y", pch = 20, col = "black")
#  matlines(LCfit_NZ$years, t(mxtHat[x, ]), lty = 1, col = "black")
#  matlines(LCfit_NZ$years, t(mxtHatPU2.5[x, ]), lty = 5, col = "red")
#  matlines(LCfit_NZ$years, t(mxtHatPU97.5[x, ]), lty = 5, col = "red")
#  matlines(LCfor_NZ$years, t(mxtCentral[x, ]), lty = 4, col = "black")
#  matlines(LCsim_NZ$years, t(mxtPred2.5[x, ]), lty = 3, col = "black")
#  matlines(LCsim_NZ$years, t(mxtPred97.5[x, ]), lty = 3, col = "black")
#  matlines(LCsimPU_NZ$years, t(mxtPredPU2.5[x, ]), lty = 5, col = "red")
#  matlines(LCsimPU_NZ$years, t(mxtPredPU97.5[x, ]), lty = 5, col = "red")
#  text(1986, mxtHatPU2.5[x, "1995"], labels = c("x=40", "x=60", "x=80"))

