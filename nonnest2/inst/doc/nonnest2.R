## ----echo=FALSE, message=FALSE-------------------------------------------
## global knitr options
knitr::opts_chunk$set(fig.path='figure/nonnest2-', fig.align='center',
                      fig.show='hold', size='footnotesize', 
                      cache.path="cache/", warning=FALSE, message=FALSE)

## packages
library("lavaan")
library("nonnest2")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("nonnest2")
#  
#  # to install the development version, run
#  # library(devtools)
#  # install_github("qpsy/nonnest2")

## ------------------------------------------------------------------------
m1 <- ' visual  =~ x1 + x2 + x3 + x4
        textual =~ x4 + x5 + x6
        speed   =~ x7 + x8 + x9 '
fit1 <- cfa(m1, data=HolzingerSwineford1939)

m2 <- ' visual  =~ x1 + x2 + x3
        textual =~ x4 + x5 + x6 + x7
        speed   =~ x7 + x8 + x9 '
fit2 <- cfa(m2, data=HolzingerSwineford1939)

## ------------------------------------------------------------------------
vuongtest(fit1, fit2)

## ------------------------------------------------------------------------
icci(fit1, fit2)

## ------------------------------------------------------------------------
m1 <- ' 
       # latent variable definitions
          ind60 =~ x1 + x2 + x3
          dem60 =~ y1 + a*y2 + b*y3 + c*y4
          dem65 =~ y5 + a*y6 + b*y7 + c*y8
     
       # regressions
         dem60 ~ ind60
         dem65 ~ ind60 + dem60
     
       # residual correlations
         y1 ~~ y5
         y2 ~~ y4 + y6
         y3 ~~ y7
         y4 ~~ y8
         y6 ~~ y8
'
     
fit1 <- sem(m1, data=PoliticalDemocracy)

m2 <- ' 
       # latent variable definitions
          ind60 =~ x1 + x2 + x3
          dem60 =~ y1 + a*y2 + b*y3 + c*y4
          dem65 =~ y5 + a*y6 + b*y7 + c*y8
     
       # regressions
         dem60 ~ ind60
         dem65 ~ ind60 + dem60
     
       # residual correlations
         y1 ~~ y3 + y5
         y2 ~~ y6
         y3 ~~ y7
         y4 ~~ y8
         y5 ~~ y7
'

fit2 <- sem(m2, data=PoliticalDemocracy)

## ------------------------------------------------------------------------
vuongtest(fit1, fit2)

## ------------------------------------------------------------------------
m3 <- ' 
       # latent variable definitions
          ind60 =~ x1 + x2 + x3
          dem60 =~ y1 + y2 + y3 + y4
          dem65 =~ y5 + y6 + y7 + y8
     
       # regressions
         dem60 ~ ind60
         dem65 ~ ind60 + dem60
     
       # residual correlations
         y1 ~~ y5
         y2 ~~ y4 + y6
         y3 ~~ y7
         y4 ~~ y8
         y6 ~~ y8
'
     
fit3 <- sem(m3, data=PoliticalDemocracy)

## ------------------------------------------------------------------------
vuongtest(fit1, fit3, nested=TRUE)

## ------------------------------------------------------------------------
anova(fit1, fit3)

