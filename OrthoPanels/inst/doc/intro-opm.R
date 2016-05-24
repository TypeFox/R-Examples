## ----knitr_setup, echo = FALSE-------------------------------------------
knitr::opts_chunk$set(tidy = TRUE, warning = FALSE)
library(OrthoPanels)

## ------------------------------------------------------------------------
head(BES_panel)

## ------------------------------------------------------------------------
BES.opm <- opm(Approve~Econ+Clegg+Brown+Cameron-NHS+Terror+PID-Tax,
data=BES_panel, index = c('n', 't'), n.samp=1000, add.time.indicators=TRUE)

## ------------------------------------------------------------------------
str(BES.opm)

## ------------------------------------------------------------------------
summary(BES.opm)

## ------------------------------------------------------------------------
confint(BES.opm, level=0.9)

## ------------------------------------------------------------------------
caterplot(BES.opm, main = "BES 2010 opm parameter estimates")
abline(v=0)

## ------------------------------------------------------------------------
plot(BES.opm, 'rho')

## ------------------------------------------------------------------------
quantile(BES.opm$samples$beta[,1]/(1-BES.opm$samples$rho), probs=c(0.025, 0.5, 0.975))

