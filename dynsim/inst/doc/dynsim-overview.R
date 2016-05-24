## ----include=FALSE--------------------------------------------------
library(dynsim)

options(prompt = "R> ", continue = "+  ",
        width = 70, useFancyQuotes = FALSE)

knitr::opts_chunk$set(fig.align='center', prompt=TRUE,
                highlight=FALSE, background="white")

## ----tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE------------
data(grunfeld, package = "dynsim")

## ----tidy=FALSE, message=FALSE--------------------------------------
library(DataCombine)

grunfeld <- slide(grunfeld, Var = "invest", GroupVar = "company", 
                  TimeVar = "year", NewVar = "InvestLag")

## ----tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE------------
M1 <- lm(invest ~ InvestLag + mvalue + kstock, data = grunfeld)

## ----tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE------------
attach(grunfeld)
Scen1 <- data.frame(InvestLag = mean(InvestLag, na.rm = TRUE),
                    mvalue = quantile(mvalue, 0.95),
                    kstock = quantile(kstock, 0.95))
Scen2 <- data.frame(InvestLag = mean(InvestLag, na.rm = TRUE),
                    mvalue = mean(mvalue),
                    kstock = mean(kstock))
Scen3 <- data.frame(InvestLag = mean(InvestLag, na.rm = TRUE),
                    mvalue = quantile(mvalue, 0.05),
                    kstock = quantile(kstock, 0.05))
detach(grunfeld)

ScenComb <- list(Scen1, Scen2, Scen3)

## ----tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE------------
Sim1 <- dynsim(obj = M1, ldv = "InvestLag", scen = ScenComb, n = 20)

## ----tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE------------
# Keep only the mvalue for the first company for the first 15 years
grunfeldsub <- subset(grunfeld, company == 1)
grunfeldshock <- grunfeldsub[1:15, "mvalue"]

# Create data frame for the shock argument
grunfeldshock <- data.frame(times = 1:15, mvalue = grunfeldshock)

## ----tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE------------
Sim2 <- dynsim(obj = M1, ldv = "InvestLag", scen = ScenComb, n = 15,
               shocks = grunfeldshock)

## ----tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE------------
M2 <- lm(invest ~ InvestLag + mvalue*kstock, data = grunfeld)

## ----tidy=FALSE, echo=TRUE, message=FALSE, warning=FALSE------------
Sim3 <- dynsim(obj = M2, ldv = "InvestLag", scen = ScenComb, n = 15,
               shocks = grunfeldshock)

## ----SimDefault, tidy=FALSE, eval=FALSE-----------------------------
#  dynsimGG(Sim1)

## ----Sim1Labels,tidy=FALSE, eval=FALSE------------------------------
#  Labels <- c("95th Percentile", "Mean", "5th Percentile")
#  
#  dynsimGG(Sim1, leg.name = "Scenarios", leg.labels = Labels, color = "OrRd",
#           ylab = "Predicted Real Gross Investment\n")

## ----Sim1Labels-eval,tidy=FALSE, echo=FALSE, message=FALSE, fig.height=3.5, fig.width=5.5----
Labels <- c("95th Percentile", "Mean", "5th Percentile")

dynsimGG(Sim1, leg.name = "Scenarios", leg.labels = Labels, color = "OrRd",
         ylab = "Predicted Real Gross Investment\n")

## ----Sim2Shock, tidy=FALSE, eval=FALSE------------------------------
#  dynsimGG(Sim2, leg.name = "Scenarios", leg.labels = Labels, color = "OrRd",
#         ylab = "Predicted Real Gross Investment\n", shockplot.var = "mvalue",
#         shockplot.ylab = "Firm Value")

## ----Sim2Shock-eval, tidy=FALSE, echo=FALSE, fig.height=5, fig.width=6----
dynsimGG(Sim2, leg.name = "Scenarios", leg.labels = Labels, color = "OrRd",
       ylab = "Predicted Real Gross Investment\n", shockplot.var = "mvalue",
       shockplot.ylab = "Firm Value", shockplot.heights = c(8, 4))

