## ----, echo = FALSE------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----, echo=TRUE---------------------------------------------------------
  suppressPackageStartupMessages(library(easyVerification))
  ls(pos="package:easyVerification")

## ----, eval=FALSE--------------------------------------------------------
#  devtools::install_github("MeteoSwiss/easyVerification")

## ----, echo=TRUE---------------------------------------------------------

signal <- outer(sin(seq(0,3*pi, length=100)) + seq(-2,2,length=100), rnorm(15), '+')
fcst <- array(rnorm(100*15*51), c(100, 15, 51)) + c(signal)
obs <- array(rnorm(100*15), c(100, 15)) + c(signal)

## ----, echo=TRUE---------------------------------------------------------
f.crps <- veriApply("FairCrps", fcst=fcst, obs=obs)

## ----, echo=TRUE---------------------------------------------------------
f.rps <- veriApply("FairRps", fcst=fcst, obs=obs, prob=c(1/3,2/3))
f.rps2 <- veriApply("FairRps", fcst=fcst, obs=obs, threshold=1)

## ----, echo=TRUE---------------------------------------------------------
f.crpss <- veriApply("FairCrpss", fcst=fcst, obs=obs)
mode(f.crpss)
names(f.crpss)
range(f.crpss$crpss)

## ----, echo=TRUE---------------------------------------------------------
fcst.ref <- array(cbind(0, obs[,-ncol(obs)]), c(dim(obs), 1))
f.crpss2 <- veriApply("FairCrpss", fcst=fcst, obs=obs, fcst.ref=fcst.ref)
par(mar=c(5,5,1,1))
plot(f.crpss$crpss, f.crpss2$crpss, asp=1,
     xlab='CRPSS against climatology', ylab='CRPSS against persistence')
grid()
abline(c(0,1))


## ----, echo=TRUE---------------------------------------------------------
fcst <- array(rnorm(prod(4:8)), 4:8)
obs <- array(rnorm(prod(4:7)), 4:7)
f.me <- veriApply('EnsMe', fcst=fcst, obs=obs)
dim(f.me)

## ----, echo=TRUE---------------------------------------------------------
fcst2 <- array(aperm(fcst, c(5,1,4,2,3)), c(8, 4, 7, 5*6))
obs2 <- array(aperm(obs, c(1,4,2,3)), c(4,7,5*6))
f.me2 <- veriApply('EnsMe', fcst=fcst2, obs=obs2, tdim=3, ensdim=1)
dim(f.me2)

## ----, echo=TRUE---------------------------------------------------------
range(c(f.me) - c(f.me2))

## ----, echo=TRUE---------------------------------------------------------
bestMember <- function(ens, obs){
  best <- apply(abs(ens - obs), 1, which.min)
  return(best)
}

## ----, echo=TRUE---------------------------------------------------------
f.best <- veriApply("bestMember", fcst=fcst, obs=obs)
range(f.best)

