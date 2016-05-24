### R code from vignette source 'RcppZiggurat.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
prettyVersion <- packageDescription("RcppZiggurat")$Version
prettyDate <- format(Sys.Date(), "%B %e, %Y")

library(RcppZiggurat)
## this RData file contains the pre-run results 
##
##   [1] "chires"   "norres"   "rspeed"   "stdres"   "zigspeed"
##
## from the two speed comparions (ie rspeed and zigspeed) as well the three 
## simulation results (ie norres, stdres, and chires) 
##
## code for the two speed tests is shown below, the simulations results were run as
##
##   #!/bin/bash
##   Rscript -e 'library(RcppZiggurat); options("mc.cores"=6L); resS <- RcppZiggurat:::standardTest(5e9); saveRDS(resS, "standardTest.rds")' &
##   Rscript -e 'library(RcppZiggurat); options("mc.cores"=6L); resC <- RcppZiggurat:::chisqTest(1e10); saveRDS(resC, "chisqTest.rds")' &
##   Rscript -e 'library(RcppZiggurat); options("mc.cores"=6L); resN <- RcppZiggurat:::normalTest(5e9); saveRDS(resN, "normalTest.rds")' &
##
load("RcppZiggurat.RData")


###################################################
### code chunk number 4: R_normal_RNGs (eval = FALSE)
###################################################
## library(microbenchmark)
## res <- microbenchmark({RNGkind(,"Kinderman-Ramage"); rnorm(1e6)},
##                       {RNGkind(,"Ahrens-Dieter"); rnorm(1e6)},
##                       {RNGkind(,"Box-Muller"); rnorm(1e6)},
##                       {RNGkind(,"Inversion"); rnorm(1e6)},
##                       times=100)
## levels(res$expr) <- c("KR", "AH", "BM", "Inv")
## #saveRDS(res, file="~/git/rcppziggurat/vignettes/Rspeed.rds")


###################################################
### code chunk number 5: R_normal_RNGs_Plot
###################################################
library(lattice)
rdf <- as.data.frame(rspeed)
rdf[,1] <- ordered(rdf[,1], levels=c("AH","KR","Inv","BM"), labels=c("AH","KR","Inv","BM"))
print(bwplot(time/1e6 ~ expr, rdf, ylab="Time in msec",
             main="Time for 100 times 1e6 normal draws",
             panel=function(...,box.ratio) {
                 panel.violin(..., col="lightgray", varwidth=FALSE, box.ratio=box.ratio)
             }))


###################################################
### code chunk number 6: Zigg_normal_RNGs (eval = FALSE)
###################################################
## library(microbenchmark)
## library(lattice)
## library(RcppZiggurat)
## N <- 1e6
## res <- microbenchmark(rnorm(N),          # Marsgalia and Tsang, JSS, 2000
##                       zrnorm(N),         # based on updated Burkardt implementation
##                       zrnormGSL(N),      # GSL's ziggurat by Voss
##                       zrnormQL(N),       # QuantLib variant
##                       zrnormGl(N),       # Gretl
##                       times=100)
## levels(res$expr) <- c("RInv", "Zigg", "ZiggGSL", "ZiggQL", "ZiggGretl")
## #saveRDS(res, file="~/git/rcppziggurat/vignettes/Zigspeed.rds")


###################################################
### code chunk number 7: Zigg_normal_RNGs_Plot
###################################################
library(lattice)
#res <- readRDS("~/git/rcppziggurat/vignettes/Zigspeed.rds")
zdf <- as.data.frame(zigspeed)
zdf[,1] <- ordered(zdf[,1],
                   levels=c("Zigg", "ZiggGSL", "ZiggQL", "ZiggGretl", "RInv"),
                   labels=c("Zigg", "ZiggGSL", "ZiggQL", "ZiggGretl", "RInv"))
print(bwplot(time/1e6 ~ expr, zdf, ylab="Time in msec",
             main="Time for 100 times 1e6 normal draws",
             panel=function(...,box.ratio) {
                 panel.violin(..., col="lightgray", varwidth=FALSE, box.ratio=box.ratio)
             }))


###################################################
### code chunk number 8: stdtest
###################################################
RcppZiggurat:::plotAll(stdres)


###################################################
### code chunk number 9: nortest
###################################################
RcppZiggurat:::plotAll(norres)


###################################################
### code chunk number 10: chisqtest
###################################################
RcppZiggurat:::plotChiSq(chires)


