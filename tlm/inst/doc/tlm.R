## ----echo=T, eval=F------------------------------------------------------
#  library(tlm)

## ----echo=F--------------------------------------------------------------
suppressMessages(require(tlm))

## ----eval=F--------------------------------------------------------------
#  vignette("tlm")

## ----eval=F--------------------------------------------------------------
#  help(package = "tlm")

## ----echo=T, eval=FALSE--------------------------------------------------
#  ?tlm

## ----echo=T, eval=FALSE--------------------------------------------------
#  ?MY

## ----echo=T, eval=FALSE--------------------------------------------------
#  ?effect

## ------------------------------------------------------------------------
data(imt)
dim(imt)
head(imt)
summary(imt)

## ------------------------------------------------------------------------
modimt <- tlm(y = logimt, x = age, data = imt, ypow = 0)

## ------------------------------------------------------------------------
modimt

## ------------------------------------------------------------------------
summary(modimt)

## ------------------------------------------------------------------------
MY(modimt)

## ------------------------------------------------------------------------
MY(modimt, npoints = 3)

## ------------------------------------------------------------------------
q13 <- quantile(imt$age, probs = c(1, 3)/4)
MY(modimt, x = q13)

## ------------------------------------------------------------------------
MY(modimt, x = q13, space = "transformed")

## ----eval = F------------------------------------------------------------
#  plot(modimt, type = "transformed", observed = T, xname = "Age (years)", yname = "IMT")
#  plot(modimt, observed = T, xname = "Age (years)", yname = "IMT (mm)")

## ----plotimttrans, echo = F, eval = F------------------------------------
#  plot(modimt, type = "transformed", observed = T, xname = "Age (years)", yname = "IMT")

## ----plotimt, echo = F, eval = F-----------------------------------------
#  plot(modimt, observed = T, xname = "Age (years)", yname = "IMT (mm)")

## ----echo = F, fig.width=8, fig.height=4---------------------------------
par(las = 1, mfrow = c(1, 2), mar = c(5, 5, 3, 3),  mgp = c(2.8, 0.6, 0))
plot(modimt, type = "transformed", observed = T, xname = "Age (years)", yname = "IMT")
plot(modimt, observed = T, xname = "Age (years)", yname = "IMT (mm)")

## ----plotimtdiag, eval = F-----------------------------------------------
#  plot(modimt, type = "diagnosis")

## ------------------------------------------------------------------------
effectInfo(modimt)

## ----echo = F, fig.width=8, fig.height=8---------------------------------
#par(las = 1, mfrow = c(1, 2), mar = c(5, 5, 3, 3),  mgp = c(2.8, 0.6, 0))
par(las = 1, mar = c(5, 5, 3, 3),  mgp = c(2.8, 0.6, 0))
plot(modimt, type = "diagnosis") 

## ------------------------------------------------------------------------
effect(modimt)

## ------------------------------------------------------------------------
q123 <- quantile(imt$age, probs = 1:3/4)   # quartiles
effect(modimt, x1 = q123)

## ------------------------------------------------------------------------
data(cotinine)
dim(cotinine)
head(cotinine)
summary(cotinine)

## ------------------------------------------------------------------------
modcot <- tlm(y = weight, x = logcotinine, data = cotinine, xpow = 0)

## ------------------------------------------------------------------------
summary(modcot)

## ----eval = F------------------------------------------------------------
#  plot(modcot, type = "transformed", observed = T, xname = "Cotinine", yname = "weight (kg)")
#  plot(modcot, xname = "Cotinine (ng/ml)", yname = "weight (kg)")

## ----plotcottrans, echo = F, eval = F------------------------------------
#  plot(modcot, type = "transformed", observed = T, xname = "Cotinine", yname = "weight (kg)")

## ----plotcot, echo = F, eval = F-----------------------------------------
#  plot(modcot, xname = "Cotinine (ng/ml)", yname = "weight (kg)")

## ----echo = F, fig.width=8, fig.height=4.5-------------------------------
par(las = 1, mfrow = c(1, 2), mar = c(5, 5, 4, 3),  mgp = c(2.8, 0.6, 0))
plot(modcot, type = "transformed", observed = T, xname = "Cotinine", yname = "weight (kg)")
plot(modcot, xname = "Cotinine (ng/ml)", yname = "weight (kg)")

## ------------------------------------------------------------------------
effectInfo(modcot)

## ------------------------------------------------------------------------
effect(modcot)

## ------------------------------------------------------------------------
effect(modcot, q = 10)

## ------------------------------------------------------------------------
range(cotinine$cotinine)
effect(modcot, x1 = 100, c = 200, npoints = 4)

## ------------------------------------------------------------------------
data(feld1)
dim(feld1)
head(feld1)
summary(feld1)

## ------------------------------------------------------------------------
modcat <-  tlm (y = logroom, x = logmattress, z = cat, data = feld1, ypow = 0, xpow = 0)

## ------------------------------------------------------------------------
summary(modcat)

## ----eval = F------------------------------------------------------------
#  plot(modcat, type = "transformed", observed = T, xname = "Mattress levels",
#       yname = "living room levels")
#  plot(modcat, xname = "Mattress levels", yname = "living room levels")

## ----plotcattrans, echo = F, eval = F------------------------------------
#  plot(modcat, type = "transformed", observed = T, xname = "Mattress levels", yname = "living room levels")

## ----plotcat, echo = F, eval = F-----------------------------------------
#  plot(modcat, xname = "Mattress levels", yname = "living room levels")

## ----echo = F, fig.width=8, fig.height=4.5-------------------------------
par(las = 1, mfrow = c(1, 2), mar = c(5, 5, 4, 3) + 0.1,  mgp = c(2.8, 0.6, 0))
plot(modcat, type = "transformed", observed = T, xname = "Mattress levels", yname = "living room levels")
plot(modcat, xname = "Mattress levels", yname = "living room levels")

## ------------------------------------------------------------------------
effectInfo(modcat)

## ------------------------------------------------------------------------
effect(modcat)

## ------------------------------------------------------------------------
modcat2 <-  tlm (y = logroom, x = cat, data = feld1, ypow = 0)
modcat2

## ------------------------------------------------------------------------
MY(modcat2)

## ------------------------------------------------------------------------
effect(modcat2)

## ----plotcat2, eval = F--------------------------------------------------
#  plot(modcat2, yname = "room levels")

## ----echo = F, fig.width=8, fig.height=4---------------------------------
m <- matrix(0, nrow = 2, ncol = 4)
m[, 2:3] <- 1
layout(m)
par(las = 1)
plot(modcat2, yname = "room levels")

## ------------------------------------------------------------------------
data(glucose)
dim(glucose)
head(glucose)
summary(glucose)

## ------------------------------------------------------------------------
modglucose <- tlm(y = inv2glu, x = inv12tri, data = glucose, ypow = -2, xpow = -1/2)
summary(modglucose)

## ------------------------------------------------------------------------
MY(modglucose)

## ----eval = F------------------------------------------------------------
#  plot(modglucose, type = "transformed", observed = T, xname = "Triglycerides (mg/dl)",
#       yname = "glucose (mg/dl)")
#  plot(modglucose, xname = "Triglycerides (mg/dl)", yname = "glucose (mg/dl)")

## ----plotglucotrans, echo = F, eval = F----------------------------------
#  plot(modglucose, type = "transformed", observed = T, xname = "Triglycerides (mg/dl)", yname = "glucose (mg/dl)")

## ----plotgluco, echo = F, eval = F---------------------------------------
#  plot(modglucose, xname = "Triglycerides (mg/dl)", yname = "glucose (mg/dl)")

## ----echo = F, fig.width=8, fig.height=5---------------------------------
par(las = 1, mfrow = c(1, 2), mar = c(5, 5, 4, 3) + 0.1,  mgp = c(3.6, 0.6, 0))
plot(modglucose, type = "transformed", observed = T, xname = "Triglycerides (mg/dl)", yname = "glucose (mg/dl)")
plot(modglucose, xname = "Triglycerides (mg/dl)", yname = "glucose (mg/dl)")

## ------------------------------------------------------------------------
effectInfo(modglucose)

## ----echo = F------------------------------------------------------------
# 2.5 and 97.5 percentiles of trigli (for text):
summtrigli <- quantile(glucose$trigly, probs = c(2.5, 25, 50, 75, 97.5) / 100)
trigli2.5 <- round(as.numeric(summtrigli[1]), 1)
trigli97.5 <- round(as.numeric(summtrigli[5]), 1)

## ------------------------------------------------------------------------
# Effects for an additive change in triglycerides level:
xc <- 50 * (1:5)
xc
effectXdiff <- effect(modglucose, x1 = xc)
effectXdiff

## ------------------------------------------------------------------------
# Effects for an percent change in triglycerides level:
xq <- 50 * 1.5^(0:4)
xq
effectXperc <- effect(modglucose, x1 = xq)
effectXperc

## ------------------------------------------------------------------------
modcot2 <- tlm(y = underweight, x = logcotinine, data = cotinine, xpow = 0, family = binomial)

## ------------------------------------------------------------------------
summary(modcot2)

## ------------------------------------------------------------------------
MY(modcot2)

## ----eval = F------------------------------------------------------------
#  plot(modcot2, type = "transformed", xname = "Cotinine (ng/ml) levels",
#       yname = "low birth weight")
#  plot(modcot2, xname = "Cotinine (ng/ml) levels", yname = "low birth weight")

## ----plotcot2trans, echo = F, eval = F-----------------------------------
#  plot(modcot2, type = "transformed", xname = "Cotinine (ng/ml) levels", yname = "low birth weight")

## ----plotcot2, echo = F, eval = F----------------------------------------
#  plot(modcot2, xname = "Cotinine (ng/ml) levels", yname = "low birth weight")

## ----echo = F, fig.width=8, fig.height=4.5-------------------------------
par(las = 1, mfrow = c(1, 2), mar = c(5, 5, 4, 3) + 0.1,  mgp = c(2.8, 0.6, 0))
plot(modcot2, type = "transformed", xname = "Cotinine (ng/ml) levels", yname = "low birth weight")
plot(modcot2, xname = "Cotinine (ng/ml) levels", yname = "low birth weight")

## ------------------------------------------------------------------------
effectInfo(modcot2)

## ------------------------------------------------------------------------
effect(modcot2)

## ------------------------------------------------------------------------
effect(modcot2, q = 10)

