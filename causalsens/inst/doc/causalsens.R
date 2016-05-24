## ----setup, include=FALSE, cache=FALSE--------------------------------------------------
# set global chunk options
library(knitr)
opts_chunk$set(fig.path='figure/minimal-', fig.align='center', fig.show='hold',dev='pdf', dev.args=list(family="Palatino"))
options(replace.assign=TRUE,width=90)
#render_listings()

## ----loaddata, include=TRUE, cache=FALSE, tidy=FALSE------------------------------------
library(causalsens)
data(lalonde.exp)

## ----outcomemodel, tidy=FALSE-----------------------------------------------------------
ymodel <- lm(re78 ~ treat+age + education + black + hispanic + married +
             nodegree + re74 + re75 + u74 + u75, data = lalonde.exp)
summary(ymodel)

## ----pscoremodel, tidy=FALSE------------------------------------------------------------
pmodel <- glm(treat ~ age + education + black + hispanic + married +
              nodegree + re74 + re75 + u74 + u75, data = lalonde.exp,
              family = binomial())
summary(pmodel)

## ----causalsens, tidy=FALSE-------------------------------------------------------------
alpha <- seq(-4500, 4500, by = 250)
ll.sens <- causalsens(ymodel, pmodel, ~ age + education, data = lalonde.exp,
                      alpha = alpha, confound = one.sided.att)

## ----sensplot,echo=2:3,eval=FALSE-------------------------------------------------------
#  par(mfrow=c(1,2))
#  plot(ll.sens, type = "raw", bty = "n")
#  plot(ll.sens, type = "r.squared", bty = "n")

## ----sensplot, fig.width=10, fig.height=5, fig.cap="Results from LaLonde (1986)",echo=FALSE----
par(mfrow=c(1,2))
plot(ll.sens, type = "raw", bty = "n")
plot(ll.sens, type = "r.squared", bty = "n")

